#include "Simulation.h"


MacCormack::MacCormack() : Simulation()
{
  // finds minimum spatial stepsize for timestep calculation
  double min_dim = dy;
  if (dx < dy) {
    min_dim = dx;
  }

  // quickly getting the maximum BC velocity to use for timestep calculation
  double max_boundary_speed = c;
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (fabs(boundary_v[i][j]) > max_boundary_speed) {
        max_boundary_speed = boundary_v[i][j];
      }
    }
  }

  // adding 1. here to the max speed to reproduce Borg's result
  // dt = 0.5 * min_dim / (1./mach + 1.0);
  dt = 0.5 * min_dim / max_boundary_speed;
  cout << "MacCormack timestep defined by stability criteria: " << dt << endl;

  rs = create2dArray<double>(grid_size_x, grid_size_y);
  rus = create2dArray<double>(grid_size_x, grid_size_y);
  rvs = create2dArray<double>(grid_size_x, grid_size_y);
  energy_s = create2dArray<double>(grid_size_x, grid_size_y);

  mu = create2dArray<double>(grid_size_x, grid_size_y);
  k = create2dArray<double>(grid_size_x, grid_size_y);
  tauxx = create2dArray<double>(grid_size_x, grid_size_y);
  tauyy = create2dArray<double>(grid_size_x, grid_size_y);
  tauxy = create2dArray<double>(grid_size_x, grid_size_y);
  qx = create2dArray<double>(grid_size_x, grid_size_y);
  qy = create2dArray<double>(grid_size_x, grid_size_y);

  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * cp / Pr;
    }
  }

  bool forward_diff_first = true;

}


void MacCormack::run_solver_step()
{

  /*
    Using forward differences for all predictor steps produces a slightly different
    result from using only backward differences. However the 2 are nearly the same when
    compared to using an alternating difference scheme. Supposedly an alternating scheme
    is more accurate but I currently have not been able to verify this.

    forward_diff_first should always be true in order to fully replicate Borg's result

    The alternating differences seems to be required for a stable incompressible flow around a square as described in Kundu.
  */

  if (TIMESTEP % 2 == 0) {
    forward_diff_first = true;
  } else {
    forward_diff_first = false;
  }


}
