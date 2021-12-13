#include "Simulation.h"


MacCormack::MacCormack() : Simulation() {

  dt = 0.5 * dx / (1./mach + u_lid);
  cout << "MacCormack timestep defined by stability criteria: " << dt << endl;

  rs = create2dArray<double>(grid_size_x, grid_size_y);
  us = create2dArray<double>(grid_size_x, grid_size_y);
  rus = create2dArray<double>(grid_size_x, grid_size_y);
  vs = create2dArray<double>(grid_size_x, grid_size_y);
  rvs = create2dArray<double>(grid_size_x, grid_size_y);

  a1 = dt / dx;
  a2 = dt / dy;
  a3 = dt / (dx*mach*mach);
  a4 = dt / (dy*mach*mach);
  a5 = 4. * dt / (3.*Re*dx*dx);
  a6 = dt / (Re*dy*dy);
  a7 = dt / (Re*dx*dx);
  a8 = 4.*dt / (3.*Re*dy*dy);
  a9 = dt / (12. * Re * dx*dy);
  a10 = 2.*(a5+a6);
  a11 = 2.*(a7+a8);

}

void MacCormack::free_flow_predictor(size_t i, size_t j) {
  rs[i][j] = r[i][j] - a1 * (ru[i+1][j] - ru[i][j]) - a2 * (rv[i][j+1] - rv[i][j]);

  rus[i][j] = ru[i][j] - a3 * (r[i+1][j] - r[i][j])
    - a1 * (r[i+1][j]*u[i+1][j]*u[i+1][j] - r[i][j]*u[i][j]*u[i][j])
    - a2 * (r[i][j+1]*u[i][j+1]*v[i][j+1] - r[i][j]*u[i][j]*v[i][j])
    - a10 * u[i][j]
    + a5 * (u[i+1][j] + u[i-1][j])
    + a6 * (u[i][j+1] + u[i][j-1])
    + a9 * (v[i+1][j+1] + v[i-1][j-1] - v[i+1][j-1] - v[i-1][j+1]);

  rvs[i][j] = rv[i][j] - a4 * (r[i][j+1] - r[i][j])
    - a1 * (r[i+1][j]*u[i+1][j]*v[i+1][j] - r[i][j]*u[i][j]*v[i][j])
    - a2 * (r[i][j+1]*v[i][j+1]*v[i][j+1] - r[i][j]*v[i][j]*v[i][j])
    - a11 * v[i][j]
    + a7 * (v[i+1][j] + v[i-1][j])
    + a8 * (v[i][j+1] + v[i][j-1])
    + a9 * (u[i+1][j+1] + u[i-1][j-1] - u[i+1][j-1] - u[i-1][j+1]);
}

void MacCormack::free_flow_corrector(size_t i, size_t j) {
  r[i][j] = 0.5 * ((r[i][j] + rs[i][j])
    - a1 * (rus[i][j] - rus[i-1][j])
    - a2 * (rvs[i][j] - rvs[i][j-1]));

  ru[i][j] = 0.5 * ((ru[i][j] + rus[i][j])
    - a3 * (rs[i][j] - rs[i-1][j])
    - a1 * (rs[i][j]*us[i][j]*us[i][j] - rs[i-1][j]*us[i-1][j]*us[i-1][j])
    - a2 * (rs[i][j]*us[i][j]*vs[i][j] - rs[i][j-1]*us[i][j-1]*vs[i][j-1])
    - a10 * us[i][j]
    + a5 * (us[i+1][j] + us[i-1][j])
    + a6 * (us[i][j+1] + us[i][j-1])
    + a9 * (vs[i+1][j+1] + vs[i-1][j-1] - vs[i+1][j-1] - vs[i-1][j+1]));

  rv[i][j] = 0.5 * ((rv[i][j] + rvs[i][j])
    - a4 * (rs[i][j] - rs[i][j-1])
    - a1 * (rs[i][j]*us[i][j]*vs[i][j] - rs[i-1][j]*us[i-1][j]*vs[i-1][j])
    - a2 * (rs[i][j]*vs[i][j]*vs[i][j] - rs[i][j-1]*vs[i][j-1]*vs[i][j-1])
    - a11 * vs[i][j]
    + a7 * (vs[i+1][j] + vs[i-1][j])
    + a8 * (vs[i][j+1] + vs[i][j-1])
    + a9 * (us[i+1][j+1] + us[i-1][j-1] - us[i+1][j-1] - us[i-1][j+1]));


}

void MacCormack::stationary_wall_predictor(size_t i, size_t j) {

  rus[i][j] = 0.;
  rvs[i][j] = 0.;

  // stationary left wall
  if (boundary[i-1][j] == -1)  {
    rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  }

  // stationary right wall
  else if (boundary[i+1][j] == -1) {
    rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
  }

  // stationary bottom wall
  else if (boundary[i][j-1] == -1) {
    rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
  }

  // stationary top wall
  else if (boundary[i][j+1] == -1) {
    rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
  }
}

void MacCormack::stationary_wall_corrector(size_t i, size_t j) {

  ru[i][j] = 0.;
  rv[i][j] = 0.;

  // stationary left wall
  if (boundary[i-1][j] == -1)  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  }

  // stationary right wall
  else if (boundary[i+1][j] == -1) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
  }

  // stationary bottom wall
  else if (boundary[i][j-1] == -1) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
  }

  // stationary top wall
  else if (boundary[i][j+1] == -1) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
  }

}

// The nature of the upwind difference at the lid edges causes density at the corners to increase indefinitely.
// I think the density at these corners should not differ too greatly from the neighboring gridpoints so
// doing an average instead of the original difference makes more physical sense.
void MacCormack::moving_wall_predictor(size_t i, size_t j) {

  // top moving lid
  if (boundary[i][j+1] == -1) {
    // Check if this is a corner point
    if (boundary[i-1][j] == 1 || boundary[i+1][j] == 1) {
      rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    } else {
      rs[i][j] = r[i][j] - 0.5*a1*u_lid * (r[i+1][j] - r[i-1][j]) + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    }

    rus[i][j] = u_lid * rs[i][j];
    rvs[i][j] = 0.;
  }

  // // bottom moving lid
  // else if (boundary[i][j-1] == -1) {
  //   // Check if this is a corner point
  //   if (boundary[i-1][j] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i+1][j] + r[i][j+1]);
  //   } else if (boundary[i+1][j] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i-1][j] + r[i][j+1]);
  //   } else {
  //     rs[i][j] = r[i][j] - 0.5*a1*u_lid * (r[i+1][j] - r[i-1][j]) - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
  //   }

  //   rus[i][j] = u_lid * rs[i][j];
  //   rvs[i][j] = 0.;
  // }

  // // left moving lid
  // else if (boundary[i-1][j] == -1) {
  //   // Check if this is a corner point
  //   if (boundary[i][j-1] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i+1][j] + r[i][j+1]);
  //   } else if (boundary[i][j+1] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i+1][j] + r[i][j-1]);
  //   } else {
  //     rs[i][j] = r[i][j] - 0.5*a2*u_lid * (r[i][j+1] - r[i][j-1]) - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  //   }

  //   rus[i][j] = 0.;
  //   rvs[i][j] = u_lid * rs[i][j];
  // }

  // // right moving lid
  // else if (boundary[i+1][j] == -1) {
  //   // Check if this is a corner point
  //   if (boundary[i][j-1] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i-1][j] + r[i][j+1]);
  //   } else if (boundary[i][j+1] == -1) {
  //     rs[i][j] = 1./3. * (r[i][j] + r[i-1][j] + r[i][j-1]);
  //   } else {
  //     rs[i][j] = r[i][j] - 0.5*a2*u_lid * (r[i][j+1] - r[i][j-1]) + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
  //   }

  //   rus[i][j] = 0.;
  //   rvs[i][j] = u_lid * rs[i][j];
  // }

}

void MacCormack::moving_wall_corrector(size_t i, size_t j) {

  // top moving lid
  if (boundary[i][j+1] == -1) {
    if (boundary[i-1][j] == 1 || boundary[i+1][j] == 1) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    }  else {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*u_lid * (rs[i+1][j] - rs[i-1][j]) + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));

    }
    ru[i][j] = r[i][j] * u_lid;
    rv[i][j] = 0.;
  }

  // // bottom moving lid
  // else if (boundary[i][j-1] == -1) {
  //   if (boundary[i-1][j] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i+1][j] + rs[i][j+1]);
  //   } else if (boundary[i+1][j] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i-1][j] + rs[i][j+1]);
  //   } else {
  //     r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*u_lid * (rs[i+1][j] - rs[i-1][j]) - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
  //     ru[i][j] = r[i][j] * u_lid;
  //     rv[i][j] = 0.;
  //   }
  // }

  // // left moving lid
  // else if (boundary[i-1][j] == -1) {
  //   if (boundary[i][j-1] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i+1][j] + rs[i][j+1]);
  //   } else if (boundary[i][j+1] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i+1][j] + rs[i][j-1]);
  //   } else {
  //     r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2*u_lid * (rs[i][j+1] - rs[i][j-1]) - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  //     ru[i][j] = 0.;
  //     rv[i][j] = r[i][j] * u_lid;
  //   }
  // }

  // // right moving lid
  // else if (boundary[i+1][j] == -1) {
  //   if (boundary[i][j-1] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i-1][j] + rs[i][j+1]);
  //   } else if (boundary[i][j+1] == -1) {
  //     r[i][j] = 1./3. * (rs[i][j] + rs[i-1][j] + rs[i][j-1]);
  //   } else {
  //     r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2*u_lid * (rs[i][j+1] - rs[i][j-1]) + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
  //     ru[i][j] = 0.;
  //     rv[i][j] = r[i][j] * u_lid;
  //   }
  // }
}

void MacCormack::run_solver_step() {

  // predictor step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i) {
    for (j=0; j<(grid_size_y); ++j) {

      if (boundary[i][j] == 0) {
        free_flow_predictor(i,j);
      } else if (boundary[i][j] == 1) {
        stationary_wall_predictor(i,j);
      } else if (boundary[i][j] == 2) {
        moving_wall_predictor(i,j);
      } else {
        continue;
      }

    }
  }

  // calculating starred velocities
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for ( i=0; i<grid_size_x; ++i) {
    for ( j=0; j<grid_size_y; ++j) {
      us[i][j] = rus[i][j] / rs[i][j];
      vs[i][j] = rvs[i][j] / rs[i][j];
    }
  }

  // Corrector step. Note the "+3" which shifts you to the corrector functions in the func_ptr
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<grid_size_x; ++i) {
    for (j=0; j<grid_size_y; ++j) {
      if (boundary[i][j] == 0) {
        free_flow_corrector(i,j);
      } else if (boundary[i][j] == 1) {
        stationary_wall_corrector(i,j);
      } else if (boundary[i][j] == 2) {
        moving_wall_corrector(i,j);
      } else {
        continue;
      }
    }
  }

  // Calculating new velocties
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for ( i=0; i<grid_size_x; ++i) {
    for ( j=0; j<grid_size_y; ++j) {
      u[i][j] = ru[i][j] / r[i][j];
      v[i][j] = rv[i][j] / r[i][j];
    }
  }
  // cout << u[1][65] << endl;
}
