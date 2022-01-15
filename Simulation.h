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
    unsigned int load_previous_run;


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

    void read_grid_and_init_struct(char * grid_file);
    void read_config();
    void check_residual();
    void record_speed(size_t x, size_t y);
    void save_speed_to_file();
    inline size_t s_i(size_t x, size_t y);

    void save_grid_variables();

    virtual void run_solver_step()=0;

  protected:
    size_t TIMESTEP;

    enum REGIONS 
    {
      EXTERNAL = -1,
      FREE_FLOW = 0,
      LEFT_WALL = 1,
      RIGHT_WALL = 2,
      TOP_WALL = 3,
      BOTTOM_WALL = 4,
      TOP_MOVING_LID = 5,
      STATIC = 6,
      RIGHT_OUTFLOW = 7,
      LEFT_INLET = 8,
      RIGHT_PRESSURE_OUTLET = 9,
      CORNER_POINT = 10,
      PERIODIC_Y_TOP = 11,
      PERIODIC_Y_BOTTOM = 12
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
    double ** tauxy_E,  ** tauxy_F; // xy component technically should be calculated differently for proper 2nd order accuracy

    // heat conduction terms derived from Fourier's law
    double ** qx;
    double ** qy;

    double k_constant; // constant used to derive thermal conductivity from viscosity

    double sutherland(double T);

    double ** speed;
    int ** region;

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

    // These control the forward and backwards differencing in the stencil itself. This removes the need for extra if statements.
    int leftx, rightx, upy, downy;

    void TAUXY(bool EorF, bool forward, size_t i, size_t j);
    void TAUXX(bool forward, size_t i, size_t j);
    void TAUYY(bool forward, size_t i, size_t j);
    void QX(bool forward, size_t i, size_t j);
    void QY(bool forward, size_t i, size_t j);

    void update_tau_and_q();
    void update_E_and_F();
    double calc_stencil(int component, size_t i, size_t j);
    void boundary_conditions(size_t i, size_t j);
    void BC_TOP_WALL(size_t i, size_t j);
    void BC_BOTTOM_WALL(size_t i, size_t j);
    void BC_RIGHT_WALL(size_t i, size_t j);
    void BC_LEFT_WALL(size_t i, size_t j);
    void BC_TOP_MOVING_LID(size_t i, size_t j);    
    void BC_RIGHT_OUTFLOW(size_t i, size_t j);
    void BC_LEFT_INLET(size_t i, size_t j);
    void BC_RIGHT_PRESSURE_OUTLET(size_t i, size_t j);
    void BC_CORNER_POINT(size_t i, size_t j);
    void BC_PERIODIC_Y_TOP(size_t i, size_t j);
    void BC_PERIODIC_Y_BOTTOM(size_t i, size_t j);
    void update_E_and_F_Periodic();

  public:
    MacCormack();

};


#endif
