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
  dt_dx = dt / dx;
  dt_dy = dt / dy;
  _dx = 1. / dx;
  _dy = 1. / dy;
    

  rs = create2dArray<double>(grid_size_x, grid_size_y);
  rus = create2dArray<double>(grid_size_x, grid_size_y);
  rvs = create2dArray<double>(grid_size_x, grid_size_y);
  energy_s = create2dArray<double>(grid_size_x, grid_size_y);

  mu = create2dArray<double>(grid_size_x, grid_size_y);
  k = create2dArray<double>(grid_size_x, grid_size_y);
  tauxx = create2dArray<double>(grid_size_x, grid_size_y);
  tauyy = create2dArray<double>(grid_size_x, grid_size_y);
  tauxy_E = create2dArray<double>(grid_size_x, grid_size_y);
  tauxy_F = create2dArray<double>(grid_size_x, grid_size_y);
  
  qx = create2dArray<double>(grid_size_x, grid_size_y);
  qy = create2dArray<double>(grid_size_x, grid_size_y);

  E0 = create2dArray<double>(grid_size_x, grid_size_y);
  E1 = create2dArray<double>(grid_size_x, grid_size_y);
  E2 = create2dArray<double>(grid_size_x, grid_size_y);
  E3 = create2dArray<double>(grid_size_x, grid_size_y);
  F0 = create2dArray<double>(grid_size_x, grid_size_y);
  F1 = create2dArray<double>(grid_size_x, grid_size_y);
  F2 = create2dArray<double>(grid_size_x, grid_size_y);
  F3 = create2dArray<double>(grid_size_x, grid_size_y);

  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * k_constant;

      // initializing starred variables to the initial values
      // so that they don't have to updated if they are fixed values
      rs[i][j] = r[i][j];
      rus[i][j] = ru[i][j];
      rvs[i][j] = rv[i][j];
      energy_s[i][j] = energy[i][j];
    }
  }

  bool forward_diff_first = true;

}


void MacCormack::TAUXY(bool EorF, bool forward, size_t i, size_t j)
{ 
  double dv_dx, du_dy;
  // E-version of tauxy
  if (EorF)
  {
    // forward difference of dv_dx
    if (forward)
    {
      if (region[i+1][j] != EXTERNAL) {
        dv_dx = (v[i+1][j] - v[i][j]) / dx;
      } else {
        dv_dx = (v[i][j] - v[i-1][j]) / dx;
      }
    }

    // backward difference of dv_dx
    
    else
    {
      if (region[i-1][j] != EXTERNAL) {
        dv_dx = (v[i][j] - v[i-1][j]) / dx;
      } else {
        dv_dx = (v[i+1][j] - v[i][j]) / dx;
      }
    }
    
    // central difference of du_dy
    if (region[i][j+1] == EXTERNAL) {
      du_dy = (3.*u[i][j] - 4.*u[i][j-1] + u[i][j-2]) / (2.*dy);
    } else if (region[i][j-1] == EXTERNAL) {
      du_dy = (-3.*u[i][j] + 4.*u[i][j+1] - u[i][j+2]) / (2.*dy);
    } else {
      du_dy = (u[i][j+1] - u[i][j-1]) / (2.*dy);
    }

    tauxy_E[i][j] = mu[i][j] * (du_dy + dv_dx);

  }

  // F-version of tauxy
  else 
  {
    // forward difference of du_dy
    if (forward)
    {
      if (region[i][j+1] != EXTERNAL) {
        du_dy = (u[i][j+1] - u[i][j]) / dy;
      } else {
        du_dy = (u[i][j] - u[i][j-1]) / dy;
      }
    }

    // backward difference of du_dy
    else 
    {
      if (region[i][j-1] != EXTERNAL) {
        du_dy = (u[i][j] - u[i][j-1]) / dy;
      } else {
        du_dy = (u[i][j+1] - u[i][j]) / dy;
      }
    }

    // central difference of dv_dx
    if (region[i+1][j] == EXTERNAL) {
      dv_dx = (3.*v[i][j] - 4.*v[i-1][j] + v[i-2][j]) / (2.*dx);
    } else if (region[i-1][j] == EXTERNAL) {
      dv_dx = (-3.*v[i][j] + 4.*v[i+1][j] - v[i+2][j]) / (2.*dx);
    } else {
      dv_dx = (v[i+1][j] - v[i-1][j]) / (2.*dx);
    }

    tauxy_F[i][j] = mu[i][j] * (du_dy + dv_dx);

  }
}

void MacCormack::TAUXX(bool forward, size_t i, size_t j) 
{ 
  double du_dx, dv_dy;

  // forward difference of du_dx
  if (forward)
  {
    if (region[i+1][j] != EXTERNAL) {
      du_dx = (u[i+1][j] - u[i][j]) / dx;
    } else {
      du_dx = (u[i][j] - u[i-1][j]) / dx;
    }
  }

  // backward difference of du_dx
  else 
  {
    if (region[i-1][j] != EXTERNAL) {
      du_dx = (u[i][j] - u[i-1][j]) / dx;
    } else {
      du_dx = (u[i+1][j] - u[i][j]) / dx;
    }
  }

  // central difference of dv_dy
  if (region[i][j+1] == EXTERNAL) {
    dv_dy = (3.*v[i][j] - 4.*v[i][j-1] + v[i][j-2]) / (2.*dy);
  } else if (region[i][j-1] == EXTERNAL) {
    dv_dy = (-3.*v[i][j] + 4.*v[i][j+1] - v[i][j+2]) / (2.*dy);
  } else {
    dv_dy = (v[i][j+1] - v[i][j-1]) / (2.*dy);
  }

  tauxx[i][j] = 2. / 3. * mu[i][j] * (2. * du_dx - dv_dy);

}

void MacCormack::TAUYY(bool forward, size_t i, size_t j) 
{
  double du_dx, dv_dy;

  // forward difference of dv_dy
  if (forward) 
  {
    if(region[i][j+1] != EXTERNAL) {
      dv_dy = (v[i][j+1] - v[i][j]) / dy;
    } else {
      dv_dy = (v[i][j] - v[i][j-1]) / dy;
    }
  }

  // backward difference of dv_dy
  else
  {
    if(region[i][j-1] != EXTERNAL) {
      dv_dy = (v[i][j] - v[i][j-1]) / dy;
    } else {
      dv_dy = (v[i][j+1] - v[i][j]) / dy;
    }
  }

  // central difference of du_dx
  if (region[i+1][j] == EXTERNAL) {
    du_dx = (3.*u[i][j] - 4.*u[i-1][j] + u[i-2][j]) / (2.*dx);
  } else if (region[i-1][j] == EXTERNAL) {
    du_dx = (-3.*u[i][j] + 4.*u[i+1][j] - u[i+2][j]) / (2.*dx);
  } else {
    du_dx = (u[i+1][j] - u[i-1][j]) / (2.*dx);
  }

  tauyy[i][j] = 2. / 3. * mu[i][j] * (2. * dv_dy - du_dx);

}

void MacCormack::QX(bool forward, size_t i, size_t j)
{
  double dT_dx;
  
  // forward difference of dT_dx
  if (forward)
  {
    if (region[i+1][j] != EXTERNAL){
      dT_dx = (temp[i+1][j] - temp[i][j]) / dx;
    } else {
      dT_dx = (temp[i][j] - temp[i-1][j]) / dx;
    }
  }

  // backward difference of dT_dx
  else
  {
    if (region[i-1][j] != EXTERNAL){
      dT_dx = (temp[i][j] - temp[i-1][j]) / dx;
    } else {
      dT_dx = (temp[i+1][j] - temp[i][j]) / dx;
    }
  }

  qx[i][j] = -k[i][j] * dT_dx;

}

void MacCormack::QY(bool forward, size_t i, size_t j)
{
  double dT_dy;
  
  // forward difference of dT_dy
  if (forward)
  {
    if (region[i][j+1] != EXTERNAL) {
      dT_dy = (temp[i][j+1] - temp[i][j]) / dy;
    } else {
      dT_dy = (temp[i][j] - temp[i][j-1]) / dy;
    }
  }

  // backward difference of dT_dy
  else
  {
    if (region[i][j-1] != EXTERNAL) {
      dT_dy = (temp[i][j] - temp[i][j-1]) / dy;
    } else {
      dT_dy = (temp[i][j+1] - temp[i][j]) / dy;
    }
  }

  qy[i][j] = -k[i][j] * dT_dy;

}

void MacCormack::update_tau_and_q()
{

  bool forward;
  if (predictor) {
    forward = !forward_diff_first;
  } else {
    forward = forward_diff_first;
  }

  // calculating tau and q over entire grid
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == EXTERNAL) continue;
      TAUXY(true, forward, i, j);
      TAUXY(false, forward, i, j);
      TAUXX(forward, i, j);
      TAUYY(forward, i, j);
      QX(forward, i, j);
      QY(forward, i, j);
    }
  }
}

void MacCormack::update_E_and_F() 
{
  double ** rho;

  // need to use r* in corrector step calculation
  if (predictor) {
    rho = r;
  } else {
    rho = rs;
  }
  // calculating E and F over entire grid
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == EXTERNAL) continue;
      
      E0[i][j] = rho[i][j]*u[i][j];
      E1[i][j] = rho[i][j]*u[i][j]*u[i][j] + p[i][j] - tauxx[i][j];
      E2[i][j] = rho[i][j]*u[i][j]*v[i][j] - tauxy_E[i][j];
      E3[i][j] = (energy[i][j] + p[i][j]) * u[i][j] - u[i][j]*tauxx[i][j] - v[i][j]*tauxy_E[i][j] + qx[i][j];

      F0[i][j] = rho[i][j]*v[i][j];
      F1[i][j] = rho[i][j]*u[i][j]*v[i][j] - tauxy_F[i][j];
      F2[i][j] = rho[i][j]*v[i][j]*v[i][j] + p[i][j] - tauyy[i][j];
      F3[i][j] = (energy[i][j] + p[i][j]) * v[i][j] - u[i][j]*tauxy_F[i][j] - v[i][j]*tauyy[i][j] + qy[i][j];

    }
  }
}

void MacCormack::boundary_conditions(size_t i, size_t j)
{
  
  if (region[i][j] == EXTERNAL) {
    return;
  } else if (region[i][j] == TOP_WALL) {
    BC_TOP_WALL(i, j);
  } else if (region[i][j] == BOTTOM_WALL) {
    BC_BOTTOM_WALL(i, j);
  } else if (region[i][j] == RIGHT_WALL) {
    BC_RIGHT_WALL(i, j);
  } else if (region[i][j] == LEFT_WALL) {
    BC_LEFT_WALL(i, j);
  } else if (region[i][j] == TOP_MOVING_LID) {
    BC_TOP_MOVING_LID(i, j);
  }

}

void MacCormack::BC_TOP_WALL(size_t i, size_t j)
{
  if (predictor) {
    rs[i][j] = r[i][j] - dt/dy * 0.5 * (3.*rv[i][j] - 4.*rv[i][j-1] + rv[i][j-2]);
    energy_s[i][j] = rs[i][j] * cv * temp[i][j];
  } else {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - dt/dy * 0.5 * (3.*rvs[i][j] - 4.*rvs[i][j-1] + rvs[i][j-2]));
    energy[i][j] = r[i][j] * cv * temp[i][j];
  }
}

void MacCormack::BC_BOTTOM_WALL(size_t i, size_t j)
{
  if (predictor) {
    rs[i][j] = r[i][j] + dt/dy * 0.5 * (3.*rv[i][j] - 4.*rv[i][j+1] + rv[i][j+2]);
    energy_s[i][j] = rs[i][j] * cv * temp[i][j];
  } else {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + dt/dy * 0.5 * (3.*rvs[i][j] - 4.*rvs[i][j+1] + rvs[i][j+2]));
    energy[i][j] = r[i][j] * cv * temp[i][j];
  }
}

void MacCormack::BC_RIGHT_WALL(size_t i, size_t j)
{ 
  if (predictor) {
    rs[i][j] = r[i][j] - dt/dx * 0.5 * (3.*ru[i][j] - 4.*ru[i-1][j] + ru[i-2][j]);
    energy_s[i][j] = rs[i][j] * cv * temp[i][j];
  } else {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - dt/dx * 0.5 * (3.*rus[i][j] - 4.*rus[i-1][j] + rus[i-2][j]));
    energy[i][j] = r[i][j] * cv * temp[i][j];
  }
}

void MacCormack::BC_LEFT_WALL(size_t i, size_t j)
{ 
  if (predictor) {
    rs[i][j] = r[i][j] + dt/dx * 0.5 * (3.*ru[i][j] - 4.*ru[i+1][j] + ru[i+2][j]);
    energy_s[i][j] = rs[i][j] * cv * temp[i][j];
  } else {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + dt/dx * 0.5 * (3.*rus[i][j] - 4.*rus[i+1][j] + rus[i+2][j]));
    energy[i][j] = r[i][j] * cv * temp[i][j];
  }
}

void MacCormack::BC_TOP_MOVING_LID(size_t i, size_t j)
{
  double rx_contribution, ry_contribution;
  double ** rho, ** rhov;
  if (predictor) {
    rho = r;
    rhov = rv;
  } else {
    rho = rs;
    rhov = rvs;
  }

  // corner point check
  if (region[i-1][j] != TOP_MOVING_LID) {
    rx_contribution =  0.;//dt/dx * 0.5 * u[i][j] * (3.*rho[i][j] - 4.*rho[i+1][j] + rho[i+2][j]);
  } else if (region[i+1][j] != TOP_MOVING_LID) {
    rx_contribution =  0.;//dt/dx * 0.5 * u[i][j] * (3.*rho[i][j] - 4.*rho[i-1][j] + rho[i-2][j]);
  } else {
    rx_contribution =  dt/dx * 0.5 * u[i][j] * (rho[i+1][j] - rho[i-1][j]);
  }

  ry_contribution = dt/dy * 0.5 * (3.*rhov[i][j] - 4.*rhov[i][j-1] + rhov[i][j-2]);

  if (predictor) 
  { 
    rs[i][j] = r[i][j] - rx_contribution - ry_contribution;
    energy_s[i][j] = rs[i][j] * (cv * temp[i][j] + 0.5 * u[i][j] * u[i][j]);
  } 
  
  else 
  {  
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - rx_contribution - ry_contribution);
    energy[i][j] = r[i][j] * (cv * temp[i][j] + 0.5 * u[i][j] * u[i][j]);
  }
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
  
  predictor = true;
  update_tau_and_q();
  update_E_and_F();

  // performing predictor step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == FREE_FLOW)
      { 
        if (forward_diff_first)
        {
          rs[i][j] = r[i][j] - dt/dx * (E0[i+1][j] - E0[i][j]) - dt/dy * (F0[i][j+1] - F0[i][j]);
          rus[i][j] = ru[i][j] - dt/dx * (E1[i+1][j] - E1[i][j]) - dt/dy * (F1[i][j+1] - F1[i][j]);
          rvs[i][j] = rv[i][j] - dt/dx * (E2[i+1][j] - E2[i][j]) - dt/dy * (F2[i][j+1] - F2[i][j]);
          energy_s[i][j] = energy[i][j] - dt/dx * (E3[i+1][j] - E3[i][j]) - dt/dy * (F3[i][j+1] - F3[i][j]);
        }
        else
        {
          rs[i][j] = r[i][j] - dt/dx * (E0[i][j] - E0[i-1][j]) - dt/dy * (F0[i][j] - F0[i][j-1]);
          rus[i][j] = ru[i][j] - dt/dx * (E1[i][j] - E1[i-1][j]) - dt/dy * (F1[i][j] - F1[i][j-1]);
          rvs[i][j] = rv[i][j] - dt/dx * (E2[i][j] - E2[i-1][j]) - dt/dy * (F2[i][j] - F2[i][j-1]);
          energy_s[i][j] = energy[i][j] - dt/dx * (E3[i][j] - E3[i-1][j]) - dt/dy * (F3[i][j] - F3[i][j-1]);
        }
      }
      else 
      {
        boundary_conditions(i, j);
      } 
    }
  }
  
  
  // recalculate primitive variables
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      u[i][j] = rus[i][j] / rs[i][j];
      v[i][j] = rvs[i][j] / rs[i][j];
      temp[i][j] = (energy_s[i][j] / rs[i][j] - (u[i][j] * u[i][j] + v[i][j] * v[i][j]) / 2.) / cv;
      p[i][j] = rs[i][j] * R * temp[i][j];
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * k_constant;
    }
  }

  predictor = false;
  update_tau_and_q();
  update_E_and_F();

  // performing corrector step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == FREE_FLOW)
      { 
        if (forward_diff_first)
        {
          r[i][j] = 0.5 * (r[i][j] + rs[i][j] - dt/dx * (E0[i][j] - E0[i-1][j]) - dt/dy * (F0[i][j] - F0[i][j-1]));
          ru[i][j] = 0.5 * (ru[i][j] + rus[i][j] - dt/dx * (E1[i][j] - E1[i-1][j]) - dt/dy * (F1[i][j] - F1[i][j-1]));
          rv[i][j] = 0.5 * (rv[i][j] + rvs[i][j] - dt/dx * (E2[i][j] - E2[i-1][j]) - dt/dy * (F2[i][j] - F2[i][j-1]));
          energy[i][j] = 0.5 * (energy[i][j] + energy_s[i][j] - dt/dx * (E3[i][j] - E3[i-1][j]) - dt/dy * (F3[i][j] - F3[i][j-1]));
        }
        else
        {
          r[i][j] = 0.5 * (r[i][j] + rs[i][j] - dt/dx * (E0[i+1][j] - E0[i][j]) - dt/dy * (F0[i][j+1] - F0[i][j]));
          ru[i][j] = 0.5 * (ru[i][j] + rus[i][j] - dt/dx * (E1[i+1][j] - E1[i][j]) - dt/dy * (F1[i][j+1] - F1[i][j]));
          rv[i][j] = 0.5 * (rv[i][j] + rvs[i][j] - dt/dx * (E2[i+1][j] - E2[i][j]) - dt/dy * (F2[i][j+1] - F2[i][j]));
          energy[i][j] = 0.5 * (energy[i][j] + energy_s[i][j] - dt/dx * (E3[i+1][j] - E3[i][j]) - dt/dy * (F3[i][j+1] - F3[i][j]));
        }
      }
      else 
      {
        boundary_conditions(i, j);
      } 
    }
  }

  // recalculate primitive variables
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      u[i][j] = ru[i][j] / r[i][j];
      v[i][j] = rv[i][j] / r[i][j];
      temp[i][j] = (energy[i][j] / r[i][j] - (u[i][j] * u[i][j] + v[i][j] * v[i][j]) / 2.) / cv;
      p[i][j] = r[i][j] * R * temp[i][j];
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * k_constant;
    }
  }

}
