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
grid_type * create1dArray(unsigned int size) {
  grid_type * v;
  v = (grid_type *) calloc(size, sizeof(grid_type));
  return v;
}

template <typename grid_type>
void delete_2d_Array(grid_type ** v, unsigned int sizex) {
  for(int i=0; i<sizex; i++) {
    free(v[i]);
  }
  free(v);
}

template <typename grid_type>
void delete_1d_Array(grid_type * v) {
  free(v);
}

Simulation::Simulation() {

  read_config();
  read_grid_and_init_struct();

  // printf("Framerate: %d\n", framerate);
  // printf("Grid size: %d x %d\n", grid_size_x, grid_size_y);
  // printf("Delta t: %lf\n", dt);
  // printf("Delta x: %lf\n", dx);
  // printf("Delta y: %lf\n", dy);
  // printf("Viscosity: %lf\n", mu);
  // printf("Speed of sound: %lf\n", c);
  // printf("%d threads detected\n", MAX_THREADS);
  // printf("Maximum allowed timestep by Courant stability: %lf\n", dx/c);

  if (run_graphics) {
    init_graphics();
  }
}

void Simulation::run() {

  for (TIMESTEP=0; TIMESTEP<max_run_time; ++TIMESTEP) {
    run_solver_step();

    if ((TIMESTEP+1) % 100 == 0 && tolerance != 0.) check_residual();
    if (run_graphics) {
      render();
      if (glfwWindowShouldClose(window)) break;
    }
  }

  // Handling end of sim operations
  save_speed_to_file();
  

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
  u_lid = 1.0;

  framerate = stoi(config["frame_rate"]);
  force = stod(config["force"]);
  run_graphics = stoi(config["run_graphics"]);
  tolerance = stod(config["tolerance"]);
  MAX_THREADS = omp_get_max_threads()-4;
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
  speed = create2dArray<double>(grid_size_x, grid_size_y);
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

int Simulation::init_graphics() {
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

  window = glfwCreateWindow(render_grid_size_x, render_grid_size_y, "BIG BOY", NULL, NULL);
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
  VERTEX_COUNT = 6*grid_size_x*grid_size_y;
  vertex_data = create1dArray<float>(VERTEX_COUNT);

  // only need to define x and y coordinates bc calloc
  for (int x=0; x<grid_size_x; ++x) {
    for (int y=0; y<grid_size_y; ++y) {
      vertex_data[(x*grid_size_x+y)*6 + 0] = ((float) x / (float) (grid_size_x-1)) * (1.f - -1.f) + -1.f;
      vertex_data[(x*grid_size_x+y)*6 + 1] = ((float) y / (float) (grid_size_y-1)) * (1.f - -1.f) + -1.f;
    }
  }


  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*VERTEX_COUNT, vertex_data, GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); // position attribute
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); // color attribute
  glEnableVertexAttribArray(1);

  glUseProgram(shaderProgram);
  glfwSwapInterval(framerate);

  return 0;
}

void Simulation::render() {
  if (glfwWindowShouldClose(window)) {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glfwTerminate();
    return;
  }
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
  memcpy(ptr, vertex_data, sizeof(float)*VERTEX_COUNT);
  glUnmapBuffer(GL_ARRAY_BUFFER);

  glPointSize(9);
  glBindVertexArray(VAO);
  glDrawArrays(GL_POINTS, 0, grid_size_x*grid_size_y);

  glfwSwapBuffers(window);
  glfwPollEvents();
}



void Simulation::check_residual() {
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

void Simulation::save_speed_to_file() {
  FILE * rho_fp;
  FILE * u_fp;
  FILE * v_fp;
  rho_fp = fopen("Data_Output/rho_file.dat","w");
  u_fp = fopen("Data_Output/U_file.dat","w");
  v_fp = fopen("Data_Output/V_file.dat","w");
  for (unsigned int y=0; y<grid_size_y; ++y) {
    for (unsigned int x=0; x<grid_size_x; ++x) {
      fprintf(rho_fp, "%.10lf", r[x][y]);
      fprintf(u_fp, "%.10lf", u[x][y]);
      fprintf(v_fp, "%.10lf", v[x][y]);
      if (x == (grid_size_x-1) ) {
        fprintf(u_fp, "\n" );
        fprintf(v_fp, "\n" );
        fprintf(rho_fp, "\n" );
      } else {
        fprintf(u_fp, " " );
        fprintf(v_fp, " " );
        fprintf(rho_fp, " " );
      }
    }
  }
  fclose(rho_fp);fclose(u_fp);fclose(v_fp);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
  // make sure the viewport matches the new window dimensions; note that width and
  // height will be significantly larger than specified on retina displays.
  glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window) {
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
      glfwSetWindowShouldClose(window, true);
  }
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
// inline size_t s_i(size_t x, size_t y) {
//   return (x*grid_size_y + y);
// }