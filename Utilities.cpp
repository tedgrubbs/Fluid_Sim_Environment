#include "Utilities.h"

struct Sim_Struct sim;
struct Info_Struct info_struct;

size_t TIMESTEP = 0;

using namespace std;

double MIN_RENDERABLE_SPEED = -DBL_MIN;
double MAX_RENDERABLE_SPEED = DBL_MAX;

// convert coordinate to opengl system
inline float convert_x_to_opengl(unsigned int x) {
  return x / (float) sim.grid_size_x;
}

inline float convert_y_to_opengl(unsigned int x) {
  return x / (float) sim.grid_size_y;
}

template <typename grid_type>
grid_type ** create2dArray(unsigned int sizex, unsigned int sizey) {
  grid_type ** v;
  v = (grid_type **) calloc(sizex, sizeof(grid_type *));
  for(int i=0; i<sizex; ++i) {
    v[i] = (grid_type *) calloc(sizey, sizeof(grid_type));
  }
  return v;
}

void record_speed(size_t x, size_t y) {
  sim.speed[x][y] = sqrt(pow(sim.u[x][y], 2.) + pow(sim.v[x][y], 2.));

  if (sim.u[x][y] > sim.u_max && sim.u[x][y] < MAX_RENDERABLE_SPEED) {
    sim.u_max = fabs(sim.u[x][y]);
  }

  if (sim.u[x][y] < sim.u_min && sim.u[x][y] > MIN_RENDERABLE_SPEED) {
    sim.u_min = fabs(sim.u[x][y]);
  }

  if (sim.v[x][y] > sim.v_max && sim.v[x][y] < MAX_RENDERABLE_SPEED) {
    sim.v_max = fabs(sim.v[x][y]);
  }

  if (sim.v[x][y] < sim.v_min && sim.v[x][y] > MIN_RENDERABLE_SPEED) {
    sim.v_min = fabs(sim.v[x][y]);
  }
}

void check_residual() {
  double res_sum = 0.;
  for (int i=0; i<sim.grid_size_x; ++i) {
    for (int j=0; j<sim.grid_size_y; ++j) {
      res_sum += fabs(sim.speed[i][j] - sim.residual[i][j]);
      sim.residual[i][j] = sim.speed[i][j];
    }
  }

  std::cout << "Residual at time " << TIMESTEP << " " << res_sum << std::endl;

  if (res_sum < sim.tolerance) {
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
  for (unsigned int y=0; y<sim.grid_size_y; ++y) {
    for (unsigned int x=0; x<sim.grid_size_x; ++x) {
      fprintf(vmag_fp, "%.10lf", sim.speed[x][y]);
      fprintf(u_fp, "%.10lf", sim.u[x][y]);
      fprintf(v_fp, "%.10lf", sim.v[x][y]);
      if (x == (sim.grid_size_x-1) ) {
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

void render() {

  glClearColor(0., 0., 0., 0.); // This sets the background color
  glClear(GL_COLOR_BUFFER_BIT);

  double color_val_x;
  double color_val_y;

  glBegin(GL_POINTS);

  // cout << sim.u_max << endl;

    for (unsigned int x=0; x<sim.grid_size_x; ++x) {
      for (unsigned int y=0; y<sim.grid_size_y; ++y) {

        if (sim.boundary[x][y] == 0) {
          color_val_x = (fabs(sim.u[x][y]) - sim.u_min) / (sim.u_max - sim.u_min);
          color_val_y = (fabs(sim.v[x][y]) - sim.v_min) / (sim.v_max - sim.v_min);

          if (fabs(sim.u[x][y]) < MIN_RENDERABLE_SPEED) {
            color_val_x = 0.;
          }
          if (fabs(sim.v[x][y]) < MIN_RENDERABLE_SPEED) {
            color_val_y = 0.;
          }
          if (fabs(sim.u[x][y]) > MAX_RENDERABLE_SPEED) {
            color_val_x = 1;
          }
          if (fabs(sim.v[x][y]) > MAX_RENDERABLE_SPEED) {
            color_val_y = 1;
          }

          glColor4f(color_val_x, 0., color_val_y, 1);
          glVertex2f(convert_x_to_opengl(x), convert_y_to_opengl(y));
        }
      }
    }

    // draws boundary. Have to draw this on top of other plot for it to show up properly
    glColor3ub(169, 169, 169);
    for (unsigned int x=0; x<sim.grid_size_x; ++x) {
      for (unsigned int y=0; y<sim.grid_size_y; ++y) {
        if (sim.boundary[x][y] != 0) {
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

void read_grid_and_init_struct() {

  sim.r = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.u = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.ru = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.v = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.rv = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.rs = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.us = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.rus = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.vs = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.speed = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.rvs = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);
  sim.boundary = create2dArray<int>(sim.grid_size_x, sim.grid_size_y);
  sim.residual = create2dArray<double>(sim.grid_size_x, sim.grid_size_y);


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

    sim.r[temp_xi][temp_yi] = temp_rho;
    sim.u[temp_xi][temp_yi] = temp_u;
    sim.ru[temp_xi][temp_yi] = temp_rho * temp_u;
    sim.v[temp_xi][temp_yi] = temp_v;
    sim.rv[temp_xi][temp_yi] = temp_rho * temp_v;
    sim.boundary[temp_xi][temp_yi] = temp_boundary;

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

  sim.grid_size_x = stoi(config["grid_size_x"]);
  sim.grid_size_y = stoi(config["grid_size_y"]);
  sim.dt = stod(config["dt"]);
  sim.dx = stod(config["dx"]);
  sim.dy = stod(config["dy"]);
  sim.mu = stod(config["viscosity"]);
  sim.c = stod(config["c"]);

  sim.mach = 0.1;
  sim.Re = 400.;
  sim.u_lid = 1.0;

  sim.a1 = sim.dt / sim.dx;
  sim.a2 = sim.dt / sim.dx;
  sim.a3 = sim.dt / (sim.dx*sim.mach*sim.mach);
  sim.a4 = sim.dt / (sim.dy*sim.mach*sim.mach);
  sim.a5 = 4. * sim.dt / (3.*sim.Re*sim.dx*sim.dx);
  sim.a6 = sim.dt / (sim.Re*sim.dy*sim.dy);
  sim.a7 = sim.dt / (sim.Re*sim.dx*sim.dx);
  sim.a8 = 4.*sim.dt / (3.*sim.Re*sim.dy*sim.dy);
  sim.a9 = sim.dt / (12. * sim.Re * sim.dx*sim.dy);
  sim.a10 = 2.*(sim.a5+sim.a6);
  sim.a11 = 2.*(sim.a7+sim.a8);

  info_struct.framerate = stoi(config["frame_rate"]);
  sim.force = stod(config["force"]);
  info_struct.run_graphics = stoi(config["run_graphics"]);
  sim.tolerance = stod(config["tolerance"]);
  info_struct.MAX_THREADS = omp_get_max_threads();
  info_struct.render_grid_size_x = stoi(config["render_grid_size_x"]);
  info_struct.render_grid_size_y = stoi(config["render_grid_size_y"]);
  info_struct.max_run_time = stoi(config["max_run_time"]);


  json_file.close();

}

// scaler_index
inline size_t s_i(size_t x, size_t y) {
  return (x*sim.grid_size_y + y);
}
