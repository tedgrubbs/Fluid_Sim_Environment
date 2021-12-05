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
#include <GL/freeglut.h>

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

  double * rho;
  double * u;
  double * rho_u;
  double * v;
  double * rho_v;

  double * speed;

  double * rho_s;
  double * u_s;
  double * rho_u_s;
  double * v_s;
  double * rho_v_s;

  // 2nd and 3rd components of E and F vectors as defined by MacCormack stencil
  double * E2;
  double * E3;
  double * F2;
  double * F3;

  double * E2_s;
  double * E3_s;
  double * F2_s;
  double * F3_s;

  int * boundary;

  // Max and min of each variable for easy normalization in rendering
  double rho_max;
  double u_max;
  double v_max;
  double rho_min;
  double u_min;
  double v_min;

  double force; // constant force term applied to all cells

  double * residual; // stores previous velocity field in order to calculate change between timesteps for convergence studies.
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
