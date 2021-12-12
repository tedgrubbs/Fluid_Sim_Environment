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

#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

const double MIN_RENDERABLE_SPEED = 0.;
const double MAX_RENDERABLE_SPEED = DBL_MAX;

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

class Simulation {
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
    size_t TIMESTEP;

    // Graphics functions and persistent variables
    int init_graphics();
    void render();
    unsigned int VBO, VAO;
    unsigned int VERTEX_COUNT;
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
    unsigned int grid_size_x, grid_size_y;
    double dt, dx, dy;
    double mu, c, mach, Re;

    double ** r;
    double ** u;
    double ** ru;
    double ** v;
    double ** rv;
    double ** speed;
    int ** boundary;

    double u_lid;
    double force;
    int MAX_THREADS;

  public:
    Simulation();
    void run();
};

class MacCormack : public Simulation {
  private:
    virtual void run_solver_step();

    double ** rs;
    double ** us;
    double ** rus;
    double ** vs;
    double ** rvs;

    double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11;

    // iteration variables defined here is faster
    size_t i;
    size_t j;

    // defining array of functions to handle different regions separately
    void (MacCormack:: * func_ptr[6])(size_t x, size_t y);
    void free_flow_predictor(size_t x, size_t y);
    void stationary_wall_predictor(size_t x, size_t y);
    void moving_wall_predictor(size_t x, size_t y);
    void free_flow_corrector(size_t x, size_t y);
    void stationary_wall_corrector(size_t x, size_t y);
    void moving_wall_corrector(size_t x, size_t y);

  public:
    MacCormack();

};


#endif
