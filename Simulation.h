#ifndef SIMULATION_H
#define SIMULATION_H

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

#include "include/glad/glad.h"
#include <GLFW/glfw3.h>

using namespace std;

const double MIN_RENDERABLE_SPEED = -DBL_MAX;
const double MAX_RENDERABLE_SPEED = DBL_MAX;
const double R = 287.0;
const double T0 = 288.16; // constant used for viscosity calculation via Sutherland's law

template <typename grid_type>
grid_type ** create2dArray(unsigned int sizex, unsigned int sizey);

template <typename grid_type>
grid_type * create1dArray(unsigned int size);

template <typename grid_type>
void delete_2d_Array(grid_type ** v, unsigned int sizex);

template <typename grid_type>
void delete_1d_Array(grid_type * v);

void leave_glut(unsigned char key, int xx, int yy);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

class Simulation 
{
  private:
    double rho_max;
    double u_max;
    double v_max;
    double rho_min;
    double u_min;
    double v_min;
    double ** residual; // stores previous velocity field in order to calculate change between timesteps for convergence studies.
    double tolerance;
    unsigned int framerate;
    unsigned int run_graphics;
    unsigned int render_grid_size_x;
    unsigned int render_grid_size_y;
    unsigned int max_run_time;


    // Graphics functions and persistent variables
    int init_graphics();
    void render();
    unsigned int VBO, VAO;
    unsigned int VERTEX_COUNT;
    unsigned int POINT_SIZE;
    float * vertex_data;
    GLFWwindow * window;
    const char * vertexShaderSource ="#version 330 core\n"
        "layout (location = 0) in vec3 aPos;\n"
        "layout (location = 1) in vec3 aColor;\n"
        "out vec3 ourColor;\n"
        "void main()\n"
        "{\n"
        "   gl_Position = vec4(aPos, 1.0);\n"
        "   ourColor = aColor;\n"
        "}\0";

    const char * fragmentShaderSource = "#version 330 core\n"
        "out vec4 FragColor;\n"
        "in vec3 ourColor;\n"
        "void main()\n"
        "{\n"
        "   FragColor = vec4(ourColor, 1.0f);\n"
        "}\n\0";

    void read_grid_and_init_struct();
    void read_config();
    void check_residual();
    void record_speed(size_t x, size_t y);
    void save_speed_to_file();
    inline size_t s_i(size_t x, size_t y);

    virtual void run_solver_step()=0;

  protected:
    size_t TIMESTEP;

    enum REGIONS 
    {
      EXTERNAL=-1,
      FREE_FLOW=0,
      LEFT_WALL=1,
      RIGHT_WALL=2,
      TOP_WALL=3,
      BOTTOM_WALL=4,
      
    };

    unsigned int grid_size_x, grid_size_y;
    double dt, dx, dy;
    
    double mu_global;     // dynamic viscosity at STP
    double c;             // speed of sound at STP. Currently only used to set timestep via stability criteria
    double cv;            // specific heat at constant volume
    double cp;            // specific heat at constant pressure
    double Pr;            // Prandtl Number
    double gamma;         // Ratio of specific heats

    double ** r;          // mass density
    double ** u;          // x-velocity
    double ** ru;         // x-momentum density
    double ** v;          // y-velocity
    double ** rv;         // y-momentum density
    double ** p;          // pressure
    double ** energy;     // total energy
    double ** int_energy; // internal energy
    double ** temp;       // temperature
    double ** mu;         // dynamic viscosity
    double ** k;          // thermal conductivity
    
    // relevant viscous stress tensor components
    double ** tauxx;      
    double ** tauyy;
    double ** tauxy_E,  ** tauxy_F;

    // heat conduction terms derived from Fourier's law
    double ** qx;
    double ** qy;

    double k_constant; // constant used to derive thermal conductivity from viscosity

    double sutherland(double T);

    double ** speed;
    int ** region;
    double ** boundary_v; // boundary condition velocity

    int MAX_THREADS;

  public:
    Simulation();
    void run();
};

class MacCormack : public Simulation 
{
  private:
    virtual void run_solver_step();

    bool predictor;
    double ** rs;
    double ** rus;
    double ** rvs;
    double ** energy_s;
    
    /*
      stencil vector components.
      0 = Density equation
      1 = X-momentum
      2 = Y-momentum
      3 = Energy
    */
    double ** E0, ** E1,  ** E2, ** E3;
    double ** F0, ** F1,  ** F2, ** F3;

    // constants that appear often in equations
    double dt_dx;
    double dt_dy;
    double _dx, _dy;
    
    bool forward_diff_first;

    // iteration variables defined here is faster
    size_t i;
    size_t j;

    void TAUXY(bool EorF, bool forward, size_t i, size_t j);
    void TAUXX(bool forward, size_t i, size_t j);
    void TAUYY(bool forward, size_t i, size_t j);
    void QX(bool forward, size_t i, size_t j);
    void QY(bool forward, size_t i, size_t j);

    void update_tau_and_q();
    void update_E_and_F();
    void boundary_conditions(size_t i, size_t j);
    void BC_TOP_WALL(size_t i, size_t j);
    void BC_BOTTOM_WALL(size_t i, size_t j);
    void BC_RIGHT_WALL(size_t i, size_t j);
    void BC_LEFT_WALL(size_t i, size_t j);

    // defining functions to handle different regions separately
    void free_flow_predictor(size_t x, size_t y);
    void stationary_wall_predictor(size_t x, size_t y);
    void moving_wall_predictor(size_t x, size_t y);
    void free_flow_corrector(size_t x, size_t y);
    void stationary_wall_corrector(size_t x, size_t y);
    void moving_wall_corrector(size_t x, size_t y);

    void inlet_predictor(size_t x, size_t y);
    void inlet_corrector(size_t x, size_t y);
    void outlet_predictor(size_t x, size_t y);
    void outlet_corrector(size_t x, size_t y);

    void extrapolate_out_predictor(size_t x, size_t y);
    void extrapolate_out_corrector(size_t x, size_t y);

    void static_u_predictor(size_t x, size_t y);
    void static_u_corrector(size_t x, size_t y);

    // stationary walls function for density equations derived from momentum navier stokes
    void stationary_wall_mom_predictor(size_t x, size_t y);
    void stationary_wall_mom_corrector(size_t x, size_t y);

  public:
    MacCormack();

};


#endif
