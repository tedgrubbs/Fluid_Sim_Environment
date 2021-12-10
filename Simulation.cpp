#include "Simulation.h"

template <typename grid_type>
grid_type ** create2dArray(unsigned int sizex, unsigned int sizey) {
  grid_type ** v;
  v = (grid_type **) calloc(sizex, sizeof(grid_type *));
  for(int i=0; i<sizex; ++i) {
    v[i] = (grid_type *) calloc(sizey, sizeof(grid_type));
  }
  return v;
}

template <typename grid_type>
void delete_2d_Array(grid_type ** v, unsigned int sizex) {
  for(int i=0; i<sizex; i++) {
    free(v[i]);
  }
  free(v);
}

Simulation::Simulation() {
  read_config();
  read_grid_and_init_struct();

  printf("Framerate: %d\n", framerate);
  printf("Grid size: %d x %d\n", grid_size_x, grid_size_y);
  printf("Delta t: %lf\n", dt);
  printf("Delta x: %lf\n", dx);
  printf("Delta y: %lf\n", dy);
  printf("Viscosity: %lf\n", mu);
  printf("Speed of sound: %lf\n", c);
  printf("%d threads detected\n", MAX_THREADS);
  printf("Maximum allowed timestep by Courant stability: %lf\n", dx/c);

  if (run_graphics) {
    init_graphics();
  }
}

void Simulation::run() {
  for (TIMESTEP=0; TIMESTEP<max_run_time; ++TIMESTEP) {
    update();
  }
}

void Simulation::update() {
  run_solver_step();

  if ((TIMESTEP+1) % 100 == 0 && tolerance != 0.) {
    check_residual();
  }

  if (run_graphics) {
    render();
  }

}

void Simulation::read_config() {
  map<string,string> config;

  ifstream json_file;
  json_file.open("config.json");

  string line;
  string json_var;
  string json_val;

  while (getline(json_file, line)) {
    if (line.find("\"") != string::npos) {
      json_var.assign(line.substr(line.find("\"") + 1, line.find(":") - line.find("\"") - 2));
      json_val.assign(line.substr(line.find(":") + 1, line.find(",") - line.find(":") - 1 ));
      config[json_var] = json_val;
    }
  }

  grid_size_x = stoi(config["grid_size_x"]);
  grid_size_y = stoi(config["grid_size_y"]);
  dt = stod(config["dt"]);
  dx = stod(config["dx"]);
  dy = stod(config["dy"]);
  mu = stod(config["viscosity"]);
  c = stod(config["c"]);

  mach = 0.1;
  Re = 400.;
  u_lid = 2.0;

  // a1 = dt / dx;
  // a2 = dt / dx;
  // a3 = dt / (dx*mach*mach);
  // a4 = dt / (dy*mach*mach);
  // a5 = 4. * dt / (3.*Re*dx*dx);
  // a6 = dt / (Re*dy*dy);
  // a7 = dt / (Re*dx*dx);
  // a8 = 4.*dt / (3.*Re*dy*dy);
  // a9 = dt / (12. * Re * dx*dy);
  // a10 = 2.*(a5+a6);
  // a11 = 2.*(a7+a8);

  framerate = stoi(config["frame_rate"]);
  force = stod(config["force"]);
  run_graphics = stoi(config["run_graphics"]);
  tolerance = stod(config["tolerance"]);
  MAX_THREADS = omp_get_max_threads();
  if (run_graphics) MAX_THREADS -= 2; // need to leave some cpu threads open when running graphics to make it smoother.
  render_grid_size_x = stoi(config["render_grid_size_x"]);
  render_grid_size_y = stoi(config["render_grid_size_y"]);
  max_run_time = stoi(config["max_run_time"]);

  json_file.close();
}


void Simulation::read_grid_and_init_struct() {

  r = create2dArray<double>(grid_size_x, grid_size_y);
  u = create2dArray<double>(grid_size_x, grid_size_y);
  ru = create2dArray<double>(grid_size_x, grid_size_y);
  v = create2dArray<double>(grid_size_x, grid_size_y);
  rv = create2dArray<double>(grid_size_x, grid_size_y);
  // rs = create2dArray<double>(grid_size_x, grid_size_y);
  // us = create2dArray<double>(grid_size_x, grid_size_y);
  // rus = create2dArray<double>(grid_size_x, grid_size_y);
  // vs = create2dArray<double>(grid_size_x, grid_size_y);
  speed = create2dArray<double>(grid_size_x, grid_size_y);
  // rvs = create2dArray<double>(grid_size_x, grid_size_y);
  boundary = create2dArray<int>(grid_size_x, grid_size_y);
  residual = create2dArray<double>(grid_size_x, grid_size_y);

  FILE * csv;
  csv = fopen("grid_variables.csv", "r");
  char line[4096];
  fgets(line, 4096, csv); // run once to get rid of header

  while (fgets(line, 4096, csv)) {

    int temp_xi, temp_yi, temp_boundary;
    double temp_rho, temp_u, temp_v;

    char * tok;

    tok = strtok(line, ",");
    temp_xi = atoi(tok);

    tok = strtok(NULL, ",\n");
    temp_yi = atoi(tok);

    tok = strtok(NULL, ",\n");
    temp_rho = atof(tok);

    tok = strtok(NULL, ",\n");
    temp_u = atof(tok);

    tok = strtok(NULL, ",\n");
    temp_v = atof(tok);

    tok = strtok(NULL, ",\n");
    temp_boundary = atof(tok);

    r[temp_xi][temp_yi] = temp_rho;
    u[temp_xi][temp_yi] = temp_u;
    ru[temp_xi][temp_yi] = temp_rho * temp_u;
    v[temp_xi][temp_yi] = temp_v;
    rv[temp_xi][temp_yi] = temp_rho * temp_v;
    boundary[temp_xi][temp_yi] = temp_boundary;
  }
  fclose(csv);
}

void Simulation::record_speed(size_t x, size_t y) {
  speed[x][y] = sqrt(pow(u[x][y], 2.) + pow(v[x][y], 2.));

  if (fabs(u[x][y]) > u_max && fabs(u[x][y]) < MAX_RENDERABLE_SPEED) {
    u_max = fabs(u[x][y]);
  }

  if (fabs(u[x][y]) < u_min && fabs(u[x][y]) > MIN_RENDERABLE_SPEED) {
    u_min = fabs(u[x][y]);
  }

  if (fabs(v[x][y]) > v_max && fabs(v[x][y]) < MAX_RENDERABLE_SPEED) {
    v_max = fabs(v[x][y]);
  }

  if (fabs(v[x][y]) < v_min && fabs(v[x][y]) > MIN_RENDERABLE_SPEED) {
    v_min = fabs(v[x][y]);
  }
}

void Simulation::init_graphics() {
  // GL initialization

  // double particle_x = 0.9;
  // double particle_y = 0.8;

  // glfw: initialize and configure
  // ------------------------------
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_VISIBLE, GL_TRUE);
  glfwWindowHint(GLFW_DOUBLEBUFFER, GL_TRUE);

  GLFWwindow* window = glfwCreateWindow(render_grid_size_x, render_grid_size_y, "BIG BOY", NULL, NULL);
  if (window == NULL) {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  // build and compile our shader program
  // ------------------------------------
  // vertex shader
  unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
  glCompileShader(vertexShader);
  // check for shader compile errors
  int success;
  char infoLog[512];
  glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
  if (!success)
  {
      glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
      std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
  }
  // fragment shader
  unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
  glCompileShader(fragmentShader);
  // check for shader compile errors
  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
    std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
  }
  // link shaders
  unsigned int shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);
  // check for linking errors
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
  }
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Setting up vertex buffers
  float vertex_data[6*grid_size_x*grid_size_y]; //(float *) calloc(6*grid_size_x*grid_size_y, sizeof(float));;
  // only need to define x and y coordinates
  for (int x=0; x<grid_size_x; ++x) {
    for (int y=0; y<grid_size_y; ++y) {
      vertex_data[(x*grid_size_x+y)*6 + 0] = ((float) x / (float) (grid_size_x-1)) * (1.f - -1.f) + -1.f;
      vertex_data[(x*grid_size_x+y)*6 + 1] = ((float) y / (float) (grid_size_y-1)) * (1.f - -1.f) + -1.f;
      vertex_data[(x*grid_size_x+y)*6 + 2] = 0.f;
      vertex_data[(x*grid_size_x+y)*6 + 3] = 0.f;
      vertex_data[(x*grid_size_x+y)*6 + 4] = 0.f;
      vertex_data[(x*grid_size_x+y)*6 + 5] = 0.f;
    }
  }


  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(vertex_data), vertex_data, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); // position attribute
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); // color attribute
  glEnableVertexAttribArray(1);

  glUseProgram(shaderProgram);
  glfwSwapInterval(0);


  // glDeleteVertexArrays(1, &VAO);
  // glDeleteBuffers(1, &VBO);
  // glDeleteProgram(shaderProgram);
  // glfwTerminate();
}

void Simulation::render() {
  u_max = -DBL_MAX;
  u_min = DBL_MAX;
  v_max = -DBL_MAX;
  v_min = DBL_MAX;

  // int particle_x_pos =  particle_x*grid_size_x;
  // int particle_y_pos =  particle_y*grid_size_y;
  //
  // particle_x = particle_x + dt*u[particle_x_pos][particle_y_pos];
  // particle_y = particle_y + dt*v[particle_x_pos][particle_y_pos];

  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2)
  for ( int x=0; x<grid_size_x; ++x) {
    for ( int y=0; y<grid_size_y; ++y) {
      record_speed(x,y);
    }
  }

  processInput(window);
  glClearColor(0.f, 0.f, 0.f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // updating vertex data

  float color_val_x;
  float color_val_y;
  float vertex_data[6*grid_size_x*grid_size_y];
  
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(color_val_x,color_val_y)
  for ( int x=0; x<grid_size_x; ++x) {
    for ( int y=0; y<grid_size_y; ++y) {

      color_val_x = (fabs(u[x][y]) - u_min) / (u_max - u_min);
      color_val_y = (fabs(v[x][y]) - v_min) / (v_max - v_min);

      if (fabs(u[x][y]) < MIN_RENDERABLE_SPEED) {
        color_val_x = 0.;
      }
      if (fabs(v[x][y]) < MIN_RENDERABLE_SPEED) {
        color_val_y = 0.;
      }
      if (fabs(u[x][y]) > MAX_RENDERABLE_SPEED) {
        color_val_x = 1.;
      }
      if (fabs(v[x][y]) > MAX_RENDERABLE_SPEED) {
        color_val_y = 1.;
      }

      vertex_data[(x*grid_size_x+y)*6 + 3] = color_val_x;
      vertex_data[(x*grid_size_x+y)*6 + 4] = 0.0;
      vertex_data[(x*grid_size_x+y)*6 + 5] = color_val_y;

    }
  }

  // vertex_data[(particle_x_pos*grid_size_x+particle_y_pos)*6 + 4] = 1.;

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  void *ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
  memcpy(ptr, vertex_data, sizeof(vertex_data));
  glUnmapBuffer(GL_ARRAY_BUFFER);

  glPointSize(9);
  glBindVertexArray(VAO);
  glDrawArrays(GL_POINTS, 0, grid_size_x*grid_size_y);

  glfwSwapBuffers(window);
  glfwPollEvents();
}

void Simulation::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
  // make sure the viewport matches the new window dimensions; note that width and
  // height will be significantly larger than specified on retina displays.
  glViewport(0, 0, width, height);
}

void Simulation::processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, true);
  }
}

void check_residual() {
  double res_sum = 0.;
  for (int i=0; i<grid_size_x; ++i) {
    for (int j=0; j<grid_size_y; ++j) {
      res_sum += fabs(speed[i][j] - residual[i][j]);
      residual[i][j] = speed[i][j];
    }
  }

  std::cout << "Residual at time " << TIMESTEP << " " << res_sum << std::endl;

  if (res_sum < tolerance) {
    save_speed_to_file();
    std::cout << "Convergence found. Exiting.\n";
    exit(0);
  }
}

void save_speed_to_file() {
  FILE * vmag_fp;
  FILE * u_fp;
  FILE * v_fp;
  vmag_fp = fopen("Data_Output/V_mag_file.csv","w");
  u_fp = fopen("Data_Output/U_file.csv","w");
  v_fp = fopen("Data_Output/V_file.csv","w");
  for (unsigned int y=0; y<grid_size_y; ++y) {
    for (unsigned int x=0; x<grid_size_x; ++x) {
      fprintf(vmag_fp, "%.10lf", speed[x][y]);
      fprintf(u_fp, "%.10lf", u[x][y]);
      fprintf(v_fp, "%.10lf", v[x][y]);
      if (x == (grid_size_x-1) ) {
        fprintf(u_fp, "\n" );
        fprintf(v_fp, "\n" );
        fprintf(vmag_fp, "\n" );
      } else {
        fprintf(u_fp, " " );
        fprintf(v_fp, " " );
        fprintf(vmag_fp, " " );
      }
    }
  }
  fclose(vmag_fp);fclose(u_fp);fclose(v_fp);
}

string remove_quotes(string s) {
  string new_s(s);
  if (new_s.find("\"") != string::npos) {
    new_s.assign(new_s.substr(new_s.find("\"") + 1, string::npos));
    new_s.assign(new_s.substr(0, new_s.find("\"")));
  }
  return new_s;
}

// scaler_index
inline size_t s_i(size_t x, size_t y) {
  return (x*grid_size_y + y);
}
