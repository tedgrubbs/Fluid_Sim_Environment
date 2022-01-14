#include "Simulation.h"

template <typename grid_type>
grid_type ** create2dArray(unsigned int sizex, unsigned int sizey) 
{
  grid_type ** v;
  v = (grid_type **) calloc(sizex, sizeof(grid_type *));

  for(int i=0; i<sizex; ++i) 
  {
    v[i] = (grid_type *) calloc(sizey, sizeof(grid_type));
  }

  return v;
}

template <typename grid_type>
grid_type * create1dArray(unsigned int size) 
{
  grid_type * v;
  v = (grid_type *) calloc(size, sizeof(grid_type));
  return v;
}

template <typename grid_type>
void delete_2d_Array(grid_type ** v, unsigned int sizex) 
{
  for(int i=0; i<sizex; i++) {
    free(v[i]);
  }
  free(v);
}

template <typename grid_type>
void delete_1d_Array(grid_type * v) 
{
  free(v);
}

Simulation::Simulation() 
{
  read_config();

  char * grid_file = "grid_variables.csv";
  if (load_previous_run) {
    grid_file = "grid_variables_state.csv";
  }
  read_grid_and_init_struct(grid_file);

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

void Simulation::run() 
{

  for (TIMESTEP=0; TIMESTEP<max_run_time; ++TIMESTEP)
  {
    // cout << TIMESTEP << endl;
    run_solver_step();

    if (run_graphics) 
    {
      render();
      if (glfwWindowShouldClose(window)) break;
    }

    if (TIMESTEP % 10000 == 0) {
      save_grid_variables();
      save_speed_to_file();
    }

  }

  // Handling end of sim operations
  save_speed_to_file();
  save_grid_variables();

}

void Simulation::save_grid_variables()
{
  FILE * outputfp;
  char * filename = "grid_variables_state.csv";
  outputfp = fopen(filename, "w");
  fprintf(outputfp, "xi,yi,rho,u,v,temperature,pressure,region\n");

  for (int x=0; x<grid_size_x; ++x) 
  {
    for (int y=0; y<grid_size_y; ++y) 
    {
      fprintf(outputfp, "%d,%d,%lf,%lf,%lf,%lf,%lf,%d\n", x, y, r[x][y], u[x][y], v[x][y], temp[x][y], p[x][y], region[x][y]);
    }
  }
  fclose(outputfp);
}

// Reads config.json to get global static simulation variables
void Simulation::read_config() 
{
  map<string,string> config;

  ifstream json_file;
  json_file.open("config.json");

  string line;
  string json_var;
  string json_val;

  while (getline(json_file, line)) 
  {
    if (line.find("\"") != string::npos) 
    {
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
  mu_global = stod(config["viscosity"]);
  c = stod(config["c"]);
  Pr = stod(config["Pr"]);
  gamma = stod(config["gamma"]);
  cv = R / (gamma-1.);
  cp = gamma * cv;
  k_constant = cp / Pr;

  framerate = stoi(config["frame_rate"]);
  run_graphics = stoi(config["run_graphics"]);
  tolerance = stod(config["tolerance"]);
  load_previous_run = stoi(config["load_previous_run"]);

  /*
    Laptop definitely runs fastest with all threads being used. More threads also helps
    the render() function update the vertex_data array faster. However this speed drops
    quickly once the cpu temperature increases.

    On my desktop which I have found that using 4 threads is the fastest option. Most likely
    because there are 4 real cores on the CPU.
  */
  MAX_THREADS = stoi(config["thread_count"]);//omp_get_max_threads();
  render_grid_size_x = stoi(config["render_grid_size_x"]);
  render_grid_size_y = stoi(config["render_grid_size_y"]);
  max_run_time = stoi(config["max_run_time"]);
  bottom = create1dArray<int>(grid_size_x);

  json_file.close();
}

// Reads in initial conditions for simulation from grid_variables.csv
void Simulation::read_grid_and_init_struct(char * grid_file) 
{ 
  r = create2dArray<double>(grid_size_x, grid_size_y);
  p = create2dArray<double>(grid_size_x, grid_size_y);
  u = create2dArray<double>(grid_size_x, grid_size_y);
  ru = create2dArray<double>(grid_size_x, grid_size_y);
  v = create2dArray<double>(grid_size_x, grid_size_y);
  rv = create2dArray<double>(grid_size_x, grid_size_y);
  temp = create2dArray<double>(grid_size_x, grid_size_y);
  energy = create2dArray<double>(grid_size_x, grid_size_y);
  int_energy = create2dArray<double>(grid_size_x, grid_size_y);
  speed = create2dArray<double>(grid_size_x, grid_size_y);
  region = create2dArray<int>(grid_size_x, grid_size_y);
  residual = create2dArray<double>(grid_size_x, grid_size_y);

  FILE * csv;
  csv = fopen(grid_file, "r");
  char line[4096];
  fgets(line, 4096, csv); // run once to get rid of header

  while (fgets(line, 4096, csv)) 
  {

    int temp_xi, temp_yi, temp_region;
    double temp_rho, temp_u, temp_v, temp_temp, temp_p;

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
    temp_temp = atof(tok);

    tok = strtok(NULL, ",\n");
    temp_p = atof(tok);

    tok = strtok(NULL, ",\n");
    temp_region = atof(tok);

    r[temp_xi][temp_yi] = temp_rho;
    p[temp_xi][temp_yi] = temp_p;
    u[temp_xi][temp_yi] = temp_u;
    ru[temp_xi][temp_yi] = temp_rho * temp_u;
    v[temp_xi][temp_yi] = temp_v;
    rv[temp_xi][temp_yi] = temp_rho * temp_v;
    temp[temp_xi][temp_yi] = temp_temp;
    int_energy[temp_xi][temp_yi] = temp_temp * cv;
    energy[temp_xi][temp_yi] = temp_rho * (int_energy[temp_xi][temp_yi] + 0.5*(temp_u*temp_u + temp_v*temp_v));
    region[temp_xi][temp_yi] = temp_region;
  }
  fclose(csv);
}

int Simulation::init_graphics() 
{
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
  if (window == NULL) 
  {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) 
  {
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
  if (!success) 
  {
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
  if (!success) 
  {
    glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
  }
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Setting up vertex buffers
  VERTEX_COUNT = 6*grid_size_x*grid_size_y;
  vertex_data = create1dArray<float>(VERTEX_COUNT);

  // only need to define x and y coordinates bc calloc
  for (int x=0; x<grid_size_x; ++x) 
  {
    for (int y=0; y<grid_size_y; ++y) 
    {
      vertex_data[(x*grid_size_y+y)*6 + 0] = ((float) x / (float) (grid_size_x-1)) * (1.f - -1.f) + -1.f;
      vertex_data[(x*grid_size_y+y)*6 + 1] = ((float) y / (float) (grid_size_y-1)) * (1.f - -1.f) + -1.f;
      if (region[x][y] == EXTERNAL) 
      {
        vertex_data[(x*grid_size_y+y)*6 + 3] = 84./255.;
        vertex_data[(x*grid_size_y+y)*6 + 4] = 84./255.;
        vertex_data[(x*grid_size_y+y)*6 + 5] = 84./255.;
      }
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

  float p_size_y = ceil((float)render_grid_size_y / (float)grid_size_y);
  float p_size_x= ceil((float)render_grid_size_x / (float)grid_size_x);
  if (p_size_x > p_size_y) {
    POINT_SIZE = p_size_x;
  } else {
    POINT_SIZE = p_size_y;
  }

  glUseProgram(shaderProgram);
  glfwSwapInterval(framerate);

  return 0;
}

void Simulation::render() 
{
  if (glfwWindowShouldClose(window)) 
  {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glfwTerminate();
    return;
  }
  u_max = -DBL_MAX;
  u_min = DBL_MAX;
  v_max = -DBL_MAX;
  v_min = DBL_MAX;
  double T_max = -DBL_MAX;
  double T_min = DBL_MAX;

  // loop indices
  int x,y;
  double absu, absv, T;
  double max_rho=0.;

  // int particle_x_pos =  particle_x*grid_size_x;
  // int particle_y_pos =  particle_y*grid_size_y;
  //
  // particle_x = particle_x + dt*u[particle_x_pos][particle_y_pos];
  // particle_y = particle_y + dt*v[particle_x_pos][particle_y_pos];

  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(x,y,absu,absv,T)
  for (x=0; x<grid_size_x; ++x) 
  {
    for (y=0; y<grid_size_y; ++y)
    {

      if (region[x][y] == EXTERNAL) {
        continue;
      }

      absu = fabs(u[x][y]);
      absv = fabs(v[x][y]);
      T = floor(temp[x][y]*1. + 0.5) / 1.;

      if (r[x][y] > max_rho) {
        max_rho = r[x][y];
      }

      if (! (r[x][y] == r[x][y]) || r[x][y] < 0.) {
        cout << "Failed from density explosion!\n";
        cout << x << " " << y << endl;
        glfwSetWindowShouldClose(window, true);
      }

      if (absu > u_max && absu < MAX_RENDERABLE_SPEED) {
        u_max = absu;
      } else if (absu < u_min && absu > MIN_RENDERABLE_SPEED) {
        u_min = absu;
      }

      if (absv > v_max && absv < MAX_RENDERABLE_SPEED) {
        v_max = absv;
      } else if (absv < v_min && absv > MIN_RENDERABLE_SPEED) {
        v_min = absv;
      }

      if (T > T_max && T < MAX_RENDERABLE_SPEED) {
        T_max = T;
      } else if (T < T_min && T > MIN_RENDERABLE_SPEED) {
        T_min = T;
      }
    }
  }
  // cout << T_min << " " << T_max << endl;
  if (max_rho > 10.) 
  {
    cout << "Failed from density explosion!\n";
    cout << max_rho << endl;
    glfwSetWindowShouldClose(window, true);
    return;
  }
  // if (T_min < 280.)
  // {
  //   cout << T_min <<  " Failed from negative temperature!\n";
  //   glfwSetWindowShouldClose(window, true);
  //   return;
  // }

  processInput(window);
  glClearColor(0.f, 0.f, 0.f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // updating vertex data

  float color_val_x;
  float color_val_y;
  float color_temp;

  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(color_val_x,color_val_y,color_temp,x,y,absu,absv,T)
  for (x=0; x<grid_size_x; ++x) 
  {
    for (y=0; y<grid_size_y; ++y) 
    {

      if (region[x][y] == EXTERNAL) {
        continue;
      }
      absu = fabs(u[x][y]);
      absv = fabs(v[x][y]);
      T = temp[x][y];

      color_val_x = (absu - u_min) / (u_max - u_min);
      color_val_y = (absv - v_min) / (v_max - v_min);
      color_temp = (T - T_min) / (T_max - T_min);

      if (absu < MIN_RENDERABLE_SPEED) {
        color_val_x = 0.;
      } else if (absu > MAX_RENDERABLE_SPEED) {
        color_val_x = 1.;
      }

      if (absv < MIN_RENDERABLE_SPEED) {
        color_val_y = 0.;
      } else if (absv > MAX_RENDERABLE_SPEED) {
        color_val_y = 1.;
      }

      if (color_temp < MIN_RENDERABLE_SPEED) {
        color_temp = 0.;
      } else if (color_temp > MAX_RENDERABLE_SPEED) {
        color_temp = 1.;
      }

      vertex_data[(x*grid_size_y+y)*6 + 3] = color_val_x;
      vertex_data[(x*grid_size_y+y)*6 + 4] = 0.0;
      vertex_data[(x*grid_size_y+y)*6 + 5] = color_val_y;

      // vertex_data[(x*grid_size_y+y)*6 + 3] = color_temp;
      // vertex_data[(x*grid_size_y+y)*6 + 4] = 0.0;
      // vertex_data[(x*grid_size_y+y)*6 + 5] = 0.0;

    }
  }

  // vertex_data[(particle_x_pos*grid_size_x+particle_y_pos)*6 + 4] = 1.;

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  void *ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
  memcpy(ptr, vertex_data, sizeof(float)*VERTEX_COUNT);
  glUnmapBuffer(GL_ARRAY_BUFFER);

  glPointSize(POINT_SIZE);
  glBindVertexArray(VAO);
  glDrawArrays(GL_POINTS, 0, grid_size_x*grid_size_y);

  glfwSwapBuffers(window);
  glfwPollEvents();
}



void Simulation::check_residual() 
{
  double res_sum = 0.;
  for (int i=0; i<grid_size_x; ++i) 
  {
    for (int j=0; j<grid_size_y; ++j) 
    {
      res_sum += fabs(speed[i][j] - residual[i][j]);
      residual[i][j] = speed[i][j];
    }
  }

  std::cout << "Residual at time " << TIMESTEP << " " << res_sum << std::endl;

  if (res_sum < tolerance) 
  {
    save_speed_to_file();
    std::cout << "Convergence found. Exiting.\n";
    exit(0);
  }
}

void Simulation::save_speed_to_file() 
{
  FILE * rho_fp;
  FILE * u_fp;
  FILE * v_fp;
  FILE * temperature_fp;
  FILE * energy_fp;
  FILE * pressure_fp;
  rho_fp = fopen("Data_Output/rho_file.dat","w");
  u_fp = fopen("Data_Output/U_file.dat","w");
  v_fp = fopen("Data_Output/V_file.dat","w");
  temperature_fp = fopen("Data_Output/temperature_file.dat","w");
  energy_fp = fopen("Data_Output/energy_file.dat","w");
  pressure_fp = fopen("Data_Output/pressure_file.dat","w");
  for (unsigned int y=0; y<grid_size_y; ++y) 
  {
    for (unsigned int x=0; x<grid_size_x; ++x) 
    {
      fprintf(rho_fp, "%.10lf", r[x][y]);
      fprintf(temperature_fp, "%.10lf", temp[x][y]);
      fprintf(energy_fp, "%.10lf", energy[x][y]);
      fprintf(pressure_fp, "%.10lf", p[x][y]);
      fprintf(u_fp, "%0.3E", u[x][y]);
      fprintf(v_fp, "%0.3E", v[x][y]);

      if (x == (grid_size_x-1)) 
      {
        fprintf(u_fp, "\n" );
        fprintf(v_fp, "\n" );
        fprintf(rho_fp, "\n" );
        fprintf(temperature_fp, "\n" );
        fprintf(energy_fp, "\n" );
        fprintf(pressure_fp, "\n" );
      } 
      else 
      {
        fprintf(u_fp, " " );
        fprintf(v_fp, " " );
        fprintf(rho_fp, " " );
        fprintf(temperature_fp, " " );
        fprintf(energy_fp, " " );
        fprintf(pressure_fp, " " );
      }
    }
  }
  fclose(rho_fp);fclose(u_fp);fclose(v_fp);fclose(temperature_fp);fclose(energy_fp);fclose(pressure_fp);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) 
{
  // make sure the viewport matches the new window dimensions; note that width and
  // height will be significantly larger than specified on retina displays.
  glViewport(0, 0, width, height);
}

void processInput(GLFWwindow *window) 
{
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) 
  {
    glfwSetWindowShouldClose(window, true);
  }
}

string remove_quotes(string s) 
{
  string new_s(s);
  if (new_s.find("\"") != string::npos) 
  {
    new_s.assign(new_s.substr(new_s.find("\"") + 1, string::npos));
    new_s.assign(new_s.substr(0, new_s.find("\"")));
  }
  return new_s;
}

double Simulation::sutherland(double T)
{
  return mu_global*pow(T/T0, 3./2.) * ((T0+110.) / (T+110.));
}

