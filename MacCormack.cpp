#include "Simulation.h"


MacCormack::MacCormack() : Simulation() {

  dt = 0.1 * dx / (1./mach);
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

  b1 = 1./3.;
  b2 = 8.*mach*mach / (9.*dx*Re);
  b3 = mach*mach / (18.*dy*Re);
  b4 = 8.*mach*mach / (9.*dy*Re);;
  b5 = mach*mach / (18.*dx*Re);;

  bool forward_diff_first;

}

void MacCormack::free_flow_predictor(size_t i, size_t j) {

  int rightx,leftx,righty,lefty;
  if (forward_diff_first) {
    rightx = i+1;
    leftx = i;
    righty = j+1;
    lefty = j;
  } else {
    rightx = i;
    leftx = i-1;
    righty = j;
    lefty = j-1;
  }

  rs[i][j] = r[i][j] - a1 * (ru[rightx][j] - ru[leftx][j]) - a2 * (rv[i][righty] - rv[i][lefty]);

  rus[i][j] = ru[i][j] - a3 * (r[rightx][j] - r[leftx][j])
    - a1 * (r[rightx][j]*u[rightx][j]*u[rightx][j] - r[leftx][j]*u[leftx][j]*u[leftx][j])
    - a2 * (r[i][righty]*u[i][righty]*v[i][righty] - r[i][lefty]*u[i][lefty]*v[i][lefty])
    - a10 * u[i][j]
    + a5 * (u[i+1][j] + u[i-1][j])
    + a6 * (u[i][j+1] + u[i][j-1])
    + a9 * (v[i+1][j+1] + v[i-1][j-1] - v[i+1][j-1] - v[i-1][j+1]);

  rvs[i][j] = rv[i][j] - a4 * (r[i][righty] - r[i][lefty])
    - a1 * (r[rightx][j]*u[rightx][j]*v[rightx][j] - r[leftx][j]*u[leftx][j]*v[leftx][j])
    - a2 * (r[i][righty]*v[i][righty]*v[i][righty] - r[i][lefty]*v[i][lefty]*v[i][lefty])
    - a11 * v[i][j]
    + a7 * (v[i+1][j] + v[i-1][j])
    + a8 * (v[i][j+1] + v[i][j-1])
    + a9 * (u[i+1][j+1] + u[i-1][j-1] - u[i+1][j-1] - u[i-1][j+1]);
}

void MacCormack::free_flow_corrector(size_t i, size_t j) {

  int rightx,leftx,righty,lefty;
  if (forward_diff_first) {
    rightx = i;
    leftx = i-1;
    righty = j;
    lefty = j-1;
  } else {
    rightx = i+1;
    leftx = i;
    righty = j+1;
    lefty = j;
  }

  r[i][j] = 0.5 * ((r[i][j] + rs[i][j])
    - a1 * (rus[rightx][j] - rus[leftx][j])
    - a2 * (rvs[i][righty] - rvs[i][lefty]));

  ru[i][j] = 0.5 * ((ru[i][j] + rus[i][j])
    - a3 * (rs[rightx][j] - rs[leftx][j])
    - a1 * (rs[rightx][j]*us[rightx][j]*us[rightx][j] - rs[leftx][j]*us[leftx][j]*us[leftx][j])
    - a2 * (rs[i][righty]*us[i][righty]*vs[i][righty] - rs[i][lefty]*us[i][lefty]*vs[i][lefty])
    - a10 * us[i][j]
    + a5 * (us[i+1][j] + us[i-1][j])
    + a6 * (us[i][j+1] + us[i][j-1])
    + a9 * (vs[i+1][j+1] + vs[i-1][j-1] - vs[i+1][j-1] - vs[i-1][j+1]));

  rv[i][j] = 0.5 * ((rv[i][j] + rvs[i][j])
    - a4 * (rs[i][righty] - rs[i][lefty])
    - a1 * (rs[rightx][j]*us[rightx][j]*vs[rightx][j] - rs[leftx][j]*us[leftx][j]*vs[leftx][j])
    - a2 * (rs[i][righty]*vs[i][righty]*vs[i][righty] - rs[i][lefty]*vs[i][lefty]*vs[i][lefty])
    - a11 * vs[i][j]
    + a7 * (vs[i+1][j] + vs[i-1][j])
    + a8 * (vs[i][j+1] + vs[i][j-1])
    + a9 * (us[i+1][j+1] + us[i-1][j-1] - us[i+1][j-1] - us[i-1][j+1]));


}

void MacCormack::stationary_wall_predictor(size_t i, size_t j) {

  rus[i][j] = 0.;
  rvs[i][j] = 0.;

  // stationary left wall
  if (region[i-1][j] == EXTERNAL)  {
    rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  }

  // stationary right wall
  else if (region[i+1][j] == EXTERNAL) {
    rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
  }

  // stationary bottom wall
  else if (region[i][j-1] == EXTERNAL) {
    rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
  }

  // stationary top wall
  else if (region[i][j+1] == EXTERNAL) {
    rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
  }
}

void MacCormack::stationary_wall_corrector(size_t i, size_t j) {

  ru[i][j] = 0.;
  rv[i][j] = 0.;

  // stationary left wall
  if (region[i-1][j] == EXTERNAL)  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  }

  // stationary right wall
  else if (region[i+1][j] == EXTERNAL) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
  }

  // stationary bottom wall
  else if (region[i][j-1] == EXTERNAL) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
  }

  // stationary top wall
  else if (region[i][j+1] == EXTERNAL) {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
  }

}

void MacCormack::moving_wall_predictor(size_t i, size_t j) {

  // top moving lid
  if (region[i][j+1] == EXTERNAL) {
    // Check if this is a corner point
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID) {
      rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    } else {
      rs[i][j] = r[i][j] - 0.5*a1*boundary_v[i][j] * (r[i+1][j] - r[i-1][j]) + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    }

    rus[i][j] = boundary_v[i][j] * rs[i][j];
    rvs[i][j] = 0.;
  }

  // bottom moving lid
  else if (region[i][j-1] == EXTERNAL) {
    // Check if this is a corner point
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID) {
      rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    } else {
      rs[i][j] = r[i][j] - 0.5*a1*boundary_v[i][j] * (r[i+1][j] - r[i-1][j]) - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    }

    rus[i][j] = boundary_v[i][j] * rs[i][j];
    rvs[i][j] = 0.;
  }



}

void MacCormack::moving_wall_corrector(size_t i, size_t j) {

  // top moving lid
  if (region[i][j+1] == EXTERNAL) {
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    }  else {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*boundary_v[i][j] * (rs[i+1][j] - rs[i-1][j]) + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));

    }
    ru[i][j] = r[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
  }

  // bottom moving lid
  else if (region[i][j-1] == EXTERNAL) {
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
    }  else {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*boundary_v[i][j] * (rs[i+1][j] - rs[i-1][j]) - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));

    }
    ru[i][j] = r[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
  }


}

void MacCormack::inlet_predictor(size_t i, size_t j) {

  // inlet on left
  if (region[i-1][j] == EXTERNAL) {
    rus[i][j] = r[i][j] * boundary_v[i][j];
    rvs[i][j] = 0.;
    rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  }

}

void MacCormack::inlet_corrector(size_t i, size_t j) {

  // inlet on left
  if (region[i-1][j] == EXTERNAL) {
    ru[i][j] = rs[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  }

}

/*
  According to some forum posts I should not be using pure outflow BCs when modeling any compressible flows or flows with varying density.
  With high reynold's number, this is a certainty and will cause the simulation to break.
  Stability is restored by maintaining a constant density at the outlet. I think this is a "pressure BC".
  Can also use velocity BC here which would make it's logic the same as the inlet. 
*/
void MacCormack::outlet_predictor(size_t i, size_t j) {

  // outlet on right
  if (region[i+1][j] == EXTERNAL) {
    rus[i][j] = 2.*ru[i-1][j] - ru[i-2][j];
    rvs[i][j] = 2.*rv[i-1][j] - rv[i-2][j];
    rs[i][j] = 1.;//2.*r[i-1][j] - r[i-2][j];
  }

}

void MacCormack::outlet_corrector(size_t i, size_t j) {

  // outlet on right
  if (region[i+1][j] == EXTERNAL) {
    ru[i][j] = 2.*rus[i-1][j] - rus[i-2][j];
    rv[i][j] = 2.*rvs[i-1][j] - rvs[i-2][j];
    r[i][j] = 1.;//2.*rs[i-1][j] - rs[i-2][j];
  }

}

void MacCormack::stationary_wall_mom_predictor(size_t i, size_t j) {

  rus[i][j] = 0.;
  rvs[i][j] = 0.;

  // if at a corner just default to normal stationary wall condition

  // left-facing wall (similar to the right wall)
  if (region[i+1][j] == EXTERNAL) {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL) {
      rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
    } else {
      rs[i][j] = b1 * (4.*r[i-1][j] - r[i-2][j])
      + b2 * (-5.*u[i-1][j] + 4.*u[i-2][j] - u[i-3][j])
      - b3 * (
        -(v[i-2][j+1] - v[i-2][j-1])
        + 4.*(v[i-1][j+1] - v[i-1][j-1])
        -3.*(v[i][j+1] - v[i][j-1])
      );
    }
  }

  // right-facing wall (similar to the left wall)
  else if(region[i-1][j] == EXTERNAL) {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL) {
      rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
    } else {
      rs[i][j] = b1 * (4.*r[i+1][j] - r[i+2][j])
      - b2 * (-5.*u[i+1][j] + 4.*u[i+2][j] - u[i+3][j])
      - b3 * (
        -(v[i+2][j+1] - v[i+2][j-1])
        + 4.*(v[i+1][j+1] - v[i+1][j-1])
        -3.*(v[i][j+1] - v[i][j-1])
      );
    }
  }

  // upward-facing wall (similar to bottom wall)
  else if (region[i][j-1] == EXTERNAL)  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL) {
      rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    } else {
      rs[i][j] = b1 * (4.*r[i][j+1] - r[i][j+2])
      - b4 * (-5.*v[i][j+1] + 4.*v[i][j+2] - v[i][j+3])
      - b5 * (
        -(u[i+1][j+2] - u[i-1][j+2])
        + 4.*(u[i+1][j+1] - u[i-1][j+1])
        -3.*(u[i+1][j] - u[i-1][j])
      );
    }
  }

  // downward facing wall (similar to top wall)
  else if (region[i][j+1] == EXTERNAL) {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL) {
      rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    } else {
      rs[i][j] = b1 * (4.*r[i][j-1] - r[i][j-2])
      + b4 * (-5.*v[i][j-1] + 4.*v[i][j-2] - v[i][j-3])
      - b5 * (
        -(u[i+1][j-2] - u[i-1][j-2])
        + 4.*(u[i+1][j-1] - u[i-1][j-1])
        -3.*(u[i+1][j] - u[i-1][j])
      );
    }
  } else {
    rs[i][j] = .25 * (r[i+1][j] + r[i-1][j] + r[i][j+1] + r[i][j-1]);
  }

}

void MacCormack::stationary_wall_mom_corrector(size_t i, size_t j) {

  ru[i][j] = 0.;
  rv[i][j] = 0.;

  // left-facing wall
  if (region[i+1][j] == EXTERNAL) {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
    } else {
      r[i][j] = b1 * (4.*rs[i-1][j] - rs[i-2][j])
      + b2 * (-5.*us[i-1][j] + 4.*us[i-2][j] - us[i-3][j])
      - b3 * (
        -(vs[i-2][j+1] - vs[i-2][j-1])
        + 4.*(vs[i-1][j+1] - vs[i-1][j-1])
        -3.*(vs[i][j+1] - vs[i][j-1])
      );
    }

  }

  // right-facing wall
  else if(region[i-1][j] == EXTERNAL) {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
    } else {
      r[i][j] = b1 * (4.*rs[i+1][j] - rs[i+2][j])
      - b2 * (-5.*us[i+1][j] + 4.*us[i+2][j] - us[i+3][j])
      - b3 * (
        -(vs[i+2][j+1] - vs[i+2][j-1])
        + 4.*(vs[i+1][j+1] - vs[i+1][j-1])
        -3.*(vs[i][j+1] - vs[i][j-1])
      );
    }

  }

  // upward-facing wall (similar to bottom wall)
  else if (region[i][j-1] == EXTERNAL)  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
    } else {
      r[i][j] = b1 * (4.*rs[i][j+1] - rs[i][j+2])
      - b4 * (-5.*vs[i][j+1] + 4.*vs[i][j+2] - vs[i][j+3])
      - b5 * (
        -(us[i+1][j+2] - us[i-1][j+2])
        + 4.*(us[i+1][j+1] - us[i-1][j+1])
        -3.*(us[i+1][j] - us[i-1][j])
      );
    }
  }

  // downward facing wall (similar to top wall)
  else if (region[i][j+1] == EXTERNAL) {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL) {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    } else {
      r[i][j] = b1 * (4.*rs[i][j-1] - rs[i][j-2])
      + b4 * (-5.*vs[i][j-1] + 4.*vs[i][j-2] - vs[i][j-3])
      - b5 * (
        -(us[i+1][j-2] - us[i-1][j-2])
        + 4.*(us[i+1][j-1] - us[i-1][j-1])
        -3.*(us[i+1][j] - us[i-1][j])
      );
    }
  } else {
    r[i][j] = .25 * (rs[i+1][j] + rs[i-1][j] + rs[i][j+1] + rs[i][j-1]);
  }

}

void MacCormack::run_solver_step() {

  /*
    Using forward differences for all predictor steps produces a slightly different
    result from using only backward differences. However the 2 are nearly the same when
    compared to using an alternating difference scheme. Supposedly an alternating scheme
    is more accurate but I currently have not been able to verify this.
  */

  if (TIMESTEP % 2 == 0) {
    forward_diff_first = true;
  } else {
    forward_diff_first = false;
  }

  // predictor step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i) {
    for (j=0; j<(grid_size_y); ++j) {

      if (region[i][j] == FREE_FLOW) {
        free_flow_predictor(i,j);
      } else if (region[i][j] == STATIONARY) {
        stationary_wall_predictor(i,j);
      } else if (region[i][j] == MOVING_LID) {
        moving_wall_predictor(i,j);
      } else if (region[i][j] == INLET) {
        inlet_predictor(i,j);
      } else if (region[i][j] == OUTLET) {
        outlet_predictor(i,j);
      } else if (region[i][j] == STATIONARY_MOM) {
        stationary_wall_mom_predictor(i,j);
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
      if (region[i][j] == FREE_FLOW) {
        free_flow_corrector(i,j);
      } else if (region[i][j] == STATIONARY) {
        stationary_wall_corrector(i,j);
      } else if (region[i][j] == MOVING_LID) {
        moving_wall_corrector(i,j);
      } else if (region[i][j] == INLET) {
        inlet_corrector(i,j);
      } else if (region[i][j] == OUTLET) {
        outlet_corrector(i,j);
      } else if (region[i][j] == STATIONARY_MOM) {
        stationary_wall_mom_corrector(i,j);
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
}
