#ifndef UTILITY
#define UTILITY

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <map>
#include <chrono>
#include <omp.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

void render();
void leave_glut(unsigned char key, int xx, int yy);
void read_grid_and_init_struct();
void read_config();
void record_speed(size_t x, size_t y);
void check_residual();
void save_speed_to_file();
inline size_t s_i(size_t x, size_t y);

// contains all relevant global variables shared across different functions
struct Sim_Struct {
  int grid_size_x;
  int grid_size_y;
  double dt;
  double dx;
  double dy;
  double mu;
  double c;
  double mach;
  double Re;
  double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11;

  double ** r;
  double ** u;
  double ** ru;
  double ** v;
  double ** rv;

  double ** speed;

  double ** rs;
  double ** us;
  double ** rus;
  double ** vs;
  double ** rvs;

  int ** boundary;
  double u_lid;

  // Max and min of each variable for easy normalization in rendering
  double rho_max;
  double u_max;
  double v_max;
  double rho_min;
  double u_min;
  double v_min;

  double force; // constant force term applied to all cells

  double ** residual; // stores previous velocity field in order to calculate change between timesteps for convergence studies.
  double tolerance;

};

struct Info_Struct {
  unsigned int framerate;
  unsigned int run_graphics;
  int MAX_THREADS;
  unsigned int render_grid_size_x;
  unsigned int render_grid_size_y;
  unsigned int max_run_time;
};

#endif
