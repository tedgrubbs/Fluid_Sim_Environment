#include "Utilities.h"
#include "Solver.h"

extern struct Sim_Struct sim_struct;
extern struct Info_Struct info_struct;

extern size_t TIMESTEP;

std::chrono::high_resolution_clock::time_point begin,end;
std::chrono::microseconds duration;

void update(int value) {



  sim_struct.u_max = -DBL_MAX;
  sim_struct.u_min = DBL_MAX;
  sim_struct.v_max = -DBL_MAX;
  sim_struct.v_min = DBL_MAX;

  run_sim_timestep();

  if ((TIMESTEP+1) % 100 == 0 && sim_struct.tolerance != 0.) {
    check_residual();
  }
  TIMESTEP += 1;

  if (info_struct.run_graphics) {
    glutPostRedisplay();
    glutTimerFunc(info_struct.framerate, update, 0);
  }

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
  printf("Grid size: %d x %d\n", sim_struct.grid_size_x, sim_struct.grid_size_y);
  printf("Delta t: %lf\n", sim_struct.dt);
  printf("Delta x: %lf\n", sim_struct.dx);
  printf("Delta y: %lf\n", sim_struct.dy);
  printf("Viscosity: %lf\n", sim_struct.mu);
  printf("Speed of sound: %lf\n", sim_struct.c);
  printf("%d threads detected\n", info_struct.MAX_THREADS);

  printf("Maximum allowed timestep by Courant stability: %lf\n", sim_struct.dx/sim_struct.c);

  if (!info_struct.run_graphics) {
    while (true) {
      update(0);
    }
  } else {
    // GL initialization
    int dummy_thicc = 0;
    glutInit(&dummy_thicc, NULL);
    glutInitDisplayMode(GLUT_SINGLE);

    //initial window size and position
    glutInitWindowSize(info_struct.render_grid_size_x, info_struct.render_grid_size_y);
    glutInitWindowPosition(1100, 200);

    //Window title and declaration of draw function
    glutCreateWindow("");
    glutDisplayFunc(render);
    glutTimerFunc(info_struct.framerate, update, 0);
    glutKeyboardFunc(leave_glut);
    glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

    float p_size_x = (float)glutGet(GLUT_WINDOW_WIDTH) / (float)sim_struct.grid_size_x ;
    float p_size_y = (float)glutGet(GLUT_WINDOW_HEIGHT) / (float)sim_struct.grid_size_y ;
    if (p_size_x > p_size_y) {
      std::cout << "Points size: " << p_size_x << std::endl;
      glPointSize(p_size_x);
    } else {
      std::cout << "Points size: " << p_size_y  << std::endl;
      glPointSize(p_size_y);
    }


    //returns you back to main() after simlation is over
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoop();
  }



  return 0;
}
