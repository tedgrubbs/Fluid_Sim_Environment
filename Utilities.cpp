#include "Utilities.h"

struct Sim_Struct sim_struct;
struct Info_Struct info_struct;

size_t TIMESTEP = 0;

using namespace std;

double MIN_RENDERABLE_SPEED = 0.;
double MAX_RENDERABLE_SPEED = 1e-6;

// convert coordinate to opengl system
inline float convert_x_to_opengl(unsigned int x) {
  return x / (float) sim_struct.grid_size_x;
}

inline float convert_y_to_opengl(unsigned int x) {
  return x / (float) sim_struct.grid_size_y;
}

void record_speed(size_t x, size_t y) {
  sim_struct.speed[s_i(x,y)] = sqrt(pow(sim_struct.u[s_i(x,y)], 2.) + pow(sim_struct.v[s_i(x,y)], 2.));

  if (sim_struct.u[s_i(x,y)] > sim_struct.u_max && sim_struct.u[s_i(x,y)] < MAX_RENDERABLE_SPEED) {
    sim_struct.u_max = fabs(sim_struct.u[s_i(x,y)]);
  }

  if (sim_struct.u[s_i(x,y)] < sim_struct.u_min && sim_struct.u[s_i(x,y)] > MIN_RENDERABLE_SPEED) {
    sim_struct.u_min = fabs(sim_struct.u[s_i(x,y)]);
  }

  if (sim_struct.v[s_i(x,y)] > sim_struct.v_max && sim_struct.v[s_i(x,y)] < MAX_RENDERABLE_SPEED) {
    sim_struct.v_max = fabs(sim_struct.v[s_i(x,y)]);
  }

  if (sim_struct.v[s_i(x,y)] < sim_struct.v_min && sim_struct.v[s_i(x,y)] > MIN_RENDERABLE_SPEED) {
    sim_struct.v_min = fabs(sim_struct.v[s_i(x,y)]);
  }
}

void check_residual() {
  double res_sum = 0.;
  for (int i=0; i<sim_struct.grid_size_x*sim_struct.grid_size_y; i++) {
    res_sum += fabs(sim_struct.speed[i] - sim_struct.residual[i]);
    sim_struct.residual[i] = sim_struct.speed[i];
  }

  std::cout << "Residual at time " << TIMESTEP << " " << res_sum << std::endl;

  if (res_sum < sim_struct.tolerance) {
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
  for (unsigned int x=0; x<sim_struct.grid_size_x; x++) {
    for (unsigned int y=0; y<sim_struct.grid_size_y; y++) {
      fprintf(vmag_fp, "%.10lf", sim_struct.speed[s_i(x,y)]);
      fprintf(u_fp, "%.10lf", sim_struct.u[s_i(x,y)]);
      fprintf(v_fp, "%.10lf", sim_struct.v[s_i(x,y)]);
      if (y == sim_struct.grid_size_y-1 && x == sim_struct.grid_size_x-1) {
        // print nothing at end
      } else {
        fprintf(vmag_fp, "," );
        fprintf(u_fp, "," );
        fprintf(v_fp, "," );
      }
    }
  }
  fclose(vmag_fp);fclose(u_fp);fclose(v_fp);
}

void render() {

  glClearColor(0., 0., 0., 0.); // This sets the background color
  glClear(GL_COLOR_BUFFER_BIT);

  double color_val_x;
  double color_val_y;

  glBegin(GL_POINTS);

  // cout << sim_struct.u_max << endl;

    for (unsigned int x=0; x<sim_struct.grid_size_x; x++) {
      for (unsigned int y=0; y<sim_struct.grid_size_y; y++) {

        if (sim_struct.boundary[s_i(x, y)] == 0) {
          color_val_x = (fabs(sim_struct.u[s_i(x, y)]) - sim_struct.u_min) / (sim_struct.u_max - sim_struct.u_min);
          color_val_y = (fabs(sim_struct.v[s_i(x, y)]) - sim_struct.v_min) / (sim_struct.v_max - sim_struct.v_min);

          if (fabs(sim_struct.u[s_i(x, y)]) < MIN_RENDERABLE_SPEED) {
            color_val_x = 0.;
          }
          if (fabs(sim_struct.v[s_i(x, y)]) < MIN_RENDERABLE_SPEED) {
            color_val_y = 0.;
          }
          if (fabs(sim_struct.u[s_i(x, y)]) > MAX_RENDERABLE_SPEED) {
            color_val_x = 1;
          }
          if (fabs(sim_struct.v[s_i(x, y)]) > MAX_RENDERABLE_SPEED) {
            color_val_y = 1;
          }

          glColor4f(color_val_x, 0., color_val_y, 1);
          glVertex2f(convert_x_to_opengl(x), convert_y_to_opengl(y));
        }
      }
    }

    // draws boundary. Have to draw this on top of other plot for it to show up properly
    glColor3ub(169, 169, 169);
    for (unsigned int x=0; x<sim_struct.grid_size_x; x++) {
      for (unsigned int y=0; y<sim_struct.grid_size_y; y++) {
        unsigned int right = (x+1) % sim_struct.grid_size_x;
        unsigned int left = (sim_struct.grid_size_x + x - 1) % sim_struct.grid_size_x;

        unsigned int up = (y+1) % sim_struct.grid_size_y;
        unsigned int down = (sim_struct.grid_size_y + y - 1) % sim_struct.grid_size_y;
        if (sim_struct.boundary[s_i(x,y)] != 0) {
          glVertex2f(convert_x_to_opengl(x), convert_y_to_opengl(y));
        }
      }
    }
  glEnd();

  glFinish();

}

void leave_glut(unsigned char key, int xx, int yy) {
  if (key == 27) {
    printf("BYEEEEEE\n");
    save_speed_to_file();
    glutLeaveMainLoop();
  }
}

double * create_double_grid() {
  return (double *) malloc(sim_struct.grid_size_x * sim_struct.grid_size_y * sizeof(double));
}

int * create_int_grid() {
  return (int *) malloc(sim_struct.grid_size_x * sim_struct.grid_size_y * sizeof(int));
}

void read_grid_and_init_struct() {

  sim_struct.rho = create_double_grid();
  sim_struct.u = create_double_grid();
  sim_struct.rho_u = create_double_grid();
  sim_struct.v = create_double_grid();
  sim_struct.rho_v = create_double_grid();
  sim_struct.rho_s = create_double_grid();
  sim_struct.u_s = create_double_grid();
  sim_struct.rho_u_s = create_double_grid();
  sim_struct.v_s = create_double_grid();
  sim_struct.speed = create_double_grid();
  sim_struct.rho_v_s = create_double_grid();
  sim_struct.E2 = create_double_grid();
  sim_struct.E3 = create_double_grid();
  sim_struct.F2 = create_double_grid();
  sim_struct.F3 = create_double_grid();
  sim_struct.E2_s = create_double_grid();
  sim_struct.E3_s = create_double_grid();
  sim_struct.F2_s = create_double_grid();
  sim_struct.F3_s = create_double_grid();
  sim_struct.boundary = create_int_grid();
  sim_struct.residual = create_double_grid();


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

    sim_struct.rho[s_i(temp_xi, temp_yi)] = temp_rho;
    sim_struct.u[s_i(temp_xi, temp_yi)] = temp_u;
    sim_struct.rho_u[s_i(temp_xi, temp_yi)] = temp_rho * temp_u;
    sim_struct.v[s_i(temp_xi, temp_yi)] = temp_v;
    sim_struct.rho_v[s_i(temp_xi, temp_yi)] = temp_rho * temp_v;
    sim_struct.boundary[s_i(temp_xi, temp_yi)] = temp_boundary;

  }
  fclose(csv);
}

string remove_quotes(string s) {
  string new_s(s);
  if (new_s.find("\"") != string::npos) {
    new_s.assign(new_s.substr(new_s.find("\"") + 1, string::npos));
    new_s.assign(new_s.substr(0, new_s.find("\"")));
  }
  return new_s;
}

void read_config() {

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

  sim_struct.grid_size_x = stoi(config["grid_size_x"]);
  sim_struct.grid_size_y = stoi(config["grid_size_y"]);
  sim_struct.dt = stod(config["dt"]);
  sim_struct.dx = stod(config["dx"]);
  sim_struct.dy = stod(config["dy"]);
  sim_struct.mu = stod(config["viscosity"]);
  sim_struct.c = stod(config["c"]);
  info_struct.framerate = stoi(config["frame_rate"]);
  sim_struct.force = stod(config["force"]);
  info_struct.run_graphics = stoi(config["run_graphics"]);
  sim_struct.tolerance = stod(config["tolerance"]);
  info_struct.MAX_THREADS = omp_get_max_threads();
  info_struct.render_grid_size_x = stoi(config["render_grid_size_x"]);
  info_struct.render_grid_size_y = stoi(config["render_grid_size_y"]);
  info_struct.max_run_time = stoi(config["max_run_time"]);


  json_file.close();

}

// scaler_index
inline size_t s_i(size_t x, size_t y) {
  return (x*sim_struct.grid_size_y + y);
}
