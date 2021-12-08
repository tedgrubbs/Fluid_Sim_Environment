#include "Utilities.h"
#include "Solver.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

const char *vertexShaderSource ="#version 330 core\n"
    "layout (location = 0) in vec3 aPos;\n"
    "layout (location = 1) in vec3 aColor;\n"
    "out vec3 ourColor;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos, 1.0);\n"
    "   ourColor = aColor;\n"
    "}\0";

const char *fragmentShaderSource = "#version 330 core\n"
    "out vec4 FragColor;\n"
    "in vec3 ourColor;\n"
    "void main()\n"
    "{\n"
    "   FragColor = vec4(ourColor, 1.0f);\n"
    "}\n\0";

extern float * vertex_data;
extern double MIN_RENDERABLE_SPEED;
extern double MAX_RENDERABLE_SPEED;

extern struct Sim_Struct sim;
extern struct Info_Struct info_struct;

extern size_t TIMESTEP;

std::chrono::high_resolution_clock::time_point begin,end;
std::chrono::microseconds duration;

void update(int value) {

  run_sim_timestep();

  if ((TIMESTEP+1) % 100 == 0 && sim.tolerance != 0.) {
    check_residual();
  }
  TIMESTEP += 1;

  if (TIMESTEP == info_struct.max_run_time) {
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << "Runtime: " << duration.count()/1000. << " milliseconds" << std::endl;
    save_speed_to_file();
    exit(0);
  }


}

int main(int argc, char const *argv[]) {

  begin = std::chrono::high_resolution_clock::now();

  read_config();
  read_grid_and_init_struct();

  printf("Framerate: %d\n", info_struct.framerate);
  printf("Grid size: %d x %d\n", sim.grid_size_x, sim.grid_size_y);
  printf("Delta t: %lf\n", sim.dt);
  printf("Delta x: %lf\n", sim.dx);
  printf("Delta y: %lf\n", sim.dy);
  printf("Viscosity: %lf\n", sim.mu);
  printf("Speed of sound: %lf\n", sim.c);
  printf("%d threads detected\n", info_struct.MAX_THREADS);

  printf("Maximum allowed timestep by Courant stability: %lf\n", sim.dx/sim.c);

  if (!info_struct.run_graphics) {
    while (true) {
      update(0);
    }
  } else {
    // GL initialization

    double particle_x = 0.9;
    double particle_y = 0.8;

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_VISIBLE, GL_TRUE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(info_struct.render_grid_size_x, info_struct.render_grid_size_y, "BIG BOY", NULL, NULL);
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
    float vertex_data[6*sim.grid_size_x*sim.grid_size_y]; //(float *) calloc(6*sim.grid_size_x*sim.grid_size_y, sizeof(float));;
    // only need to define x and y coordinates
    for (int x=0; x<sim.grid_size_x; ++x) {
      for (int y=0; y<sim.grid_size_y; ++y) {
        vertex_data[(x*sim.grid_size_x+y)*6 + 0] = ((float) x / (float) (sim.grid_size_x-1)) * (1.f - -1.f) + -1.f;
        vertex_data[(x*sim.grid_size_x+y)*6 + 1] = ((float) y / (float) (sim.grid_size_y-1)) * (1.f - -1.f) + -1.f;
        vertex_data[(x*sim.grid_size_x+y)*6 + 2] = 0.f;
        vertex_data[(x*sim.grid_size_x+y)*6 + 3] = 0.f;
        vertex_data[(x*sim.grid_size_x+y)*6 + 4] = 0.f;
        vertex_data[(x*sim.grid_size_x+y)*6 + 5] = 0.f;
      }
    }

    unsigned int VBO, VAO;
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
    while (!glfwWindowShouldClose(window)) {

      sim.u_max = -DBL_MAX;
      sim.u_min = DBL_MAX;
      sim.v_max = -DBL_MAX;
      sim.v_min = DBL_MAX;

      update(0);

      int particle_x_pos =  particle_x*sim.grid_size_x;
      int particle_y_pos =  particle_y*sim.grid_size_y;

      particle_x = particle_x + sim.dt*sim.u[particle_x_pos][particle_y_pos];
      particle_y = particle_y + sim.dt*sim.v[particle_x_pos][particle_y_pos];

      #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2)
      for ( int x=0; x<sim.grid_size_x; ++x) {
        for ( int y=0; y<sim.grid_size_y; ++y) {
          record_speed(x,y);
        }
      }

      processInput(window);
      glClearColor(0.f, 0.f, 0.f, 1.0f);
      glClear(GL_COLOR_BUFFER_BIT);

      // updating vertex data

      float color_val_x;
      float color_val_y;
      #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2) private(color_val_x,color_val_y)
      for ( int x=0; x<sim.grid_size_x; ++x) {
        for ( int y=0; y<sim.grid_size_y; ++y) {

          color_val_x = (fabs(sim.u[x][y]) - sim.u_min) / (sim.u_max - sim.u_min);
          color_val_y = (fabs(sim.v[x][y]) - sim.v_min) / (sim.v_max - sim.v_min);

          if (fabs(sim.u[x][y]) < MIN_RENDERABLE_SPEED) {
            color_val_x = 0.;
          }
          if (fabs(sim.v[x][y]) < MIN_RENDERABLE_SPEED) {
            color_val_y = 0.;
          }
          if (fabs(sim.u[x][y]) > MAX_RENDERABLE_SPEED) {
            color_val_x = 1.;
          }
          if (fabs(sim.v[x][y]) > MAX_RENDERABLE_SPEED) {
            color_val_y = 1.;
          }

          vertex_data[(x*sim.grid_size_x+y)*6 + 3] = color_val_x;
          vertex_data[(x*sim.grid_size_x+y)*6 + 4] = 0.0;
          vertex_data[(x*sim.grid_size_x+y)*6 + 5] = color_val_y;

        }
      }

      vertex_data[(particle_x_pos*sim.grid_size_x+particle_y_pos)*6 + 4] = 1.;

      glBindBuffer(GL_ARRAY_BUFFER, VBO);
      void *ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
      memcpy(ptr, vertex_data, sizeof(vertex_data));
      glUnmapBuffer(GL_ARRAY_BUFFER);

      glPointSize(9);
      glBindVertexArray(VAO);
      glDrawArrays(GL_POINTS, 0, sim.grid_size_x*sim.grid_size_y);

      glfwSwapBuffers(window);
      glfwPollEvents();

    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteProgram(shaderProgram);
    glfwTerminate();

  }

  return 0;
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
