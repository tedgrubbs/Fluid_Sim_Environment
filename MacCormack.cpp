#include "Simulation.h"

double T_STACK = 1000.;
double HEAT_RATE = 0.00712;

MacCormack::MacCormack() : Simulation()
{

  FILE * u_fp;
  u_fp = fopen("Data_Output/Probe.dat","w");
  fclose(u_fp);

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
      if (fabs(u[i][j]) > max_boundary_speed) {
        max_boundary_speed = u[i][j];
      } else if (fabs(v[i][j]) > max_boundary_speed) {
        max_boundary_speed = v[i][j];
      }
    }
  }

  // adding 1. here to the max speed to reproduce Borg's result
  // dt = 0.5 * min_dim / (1./mach + 1.0);
  dt = 0.25 * min_dim / max_boundary_speed;
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

  forward_diff_first = true;

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
  double ** rho, ** en;

  // need to use r* and E* in corrector step calculation
  if (predictor) 
  {
    rho = r;
    en = energy;
  } 
  
  else 
  {
    rho = rs;
    en = energy_s;
  }

  // calculating E and F over entire grid
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == EXTERNAL) {
        continue;
      } 
      
      E0[i][j] = rho[i][j]*u[i][j];
      E1[i][j] = rho[i][j]*u[i][j]*u[i][j] + p[i][j] - tauxx[i][j];
      E2[i][j] = rho[i][j]*u[i][j]*v[i][j] - tauxy_E[i][j];
      E3[i][j] = (en[i][j] + p[i][j]) * u[i][j] - u[i][j]*tauxx[i][j] - v[i][j]*tauxy_E[i][j] + qx[i][j];

      F0[i][j] = rho[i][j]*v[i][j];
      F1[i][j] = rho[i][j]*u[i][j]*v[i][j] - tauxy_F[i][j];
      F2[i][j] = rho[i][j]*v[i][j]*v[i][j] + p[i][j] - tauyy[i][j];
      F3[i][j] = (en[i][j] + p[i][j]) * v[i][j] - u[i][j]*tauxy_F[i][j] - v[i][j]*tauyy[i][j] + qy[i][j];

    }
  }
}

void MacCormack::update_E_and_F_Periodic()
{
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == PERIODIC_Y_TOP) 
      {
        E0[i][j] = E0[i][2];
        E1[i][j] = E1[i][2];
        E2[i][j] = E2[i][2];
        E3[i][j] = E3[i][2];

        F0[i][j] = F0[i][2];
        F1[i][j] = F1[i][2];
        F2[i][j] = F2[i][2];
        F3[i][j] = F3[i][2];
      } 
      
      else if (region[i][j] == PERIODIC_Y_BOTTOM) 
      {
        E0[i][j] = E0[i][grid_size_y-3];
        E1[i][j] = E1[i][grid_size_y-3];
        E2[i][j] = E2[i][grid_size_y-3];
        E3[i][j] = E3[i][grid_size_y-3];

        F0[i][j] = F0[i][grid_size_y-3];
        F1[i][j] = F1[i][grid_size_y-3];
        F2[i][j] = F2[i][grid_size_y-3];
        F3[i][j] = F3[i][grid_size_y-3];
      }
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
  } else if (region[i][j] == RIGHT_OUTFLOW) {
    BC_RIGHT_OUTFLOW(i, j);
  } else if (region[i][j] == LEFT_INLET) {
    BC_LEFT_INLET(i, j);
  } else if (region[i][j] == RIGHT_PRESSURE_OUTLET) {
    BC_RIGHT_PRESSURE_OUTLET(i, j);
  } else if (region[i][j] == CORNER_POINT) {
    BC_CORNER_POINT(i, j);
  } else if (region[i][j] == PERIODIC_Y_TOP) {
    BC_PERIODIC_Y_TOP(i, j);
  } else if (region[i][j] == PERIODIC_Y_BOTTOM) {
    BC_PERIODIC_Y_BOTTOM(i, j);
  }


}

/*
  Hacky way of doing periodic boundary conditions. I create a new line of grid points that take on the values of
  whatever is on the opposite side of the simulation.
*/
void MacCormack::BC_PERIODIC_Y_TOP(size_t i, size_t j)
{
  if (predictor) 
  { 
    rs[i][j] = r[i][2];
    rus[i][j] = ru[i][2];
    rvs[i][j] = rv[i][2];
    energy_s[i][j] = energy[i][2];
  }
  
  else
  {
    r[i][j] = rs[i][2];
    ru[i][j] = rus[i][2];
    rv[i][j] = rvs[i][2];
    energy[i][j] = energy_s[i][2];
  }
}

void MacCormack::BC_PERIODIC_Y_BOTTOM(size_t i, size_t j)
{
  if (predictor) 
  {
    rs[i][j] = r[i][grid_size_y-3];
    rus[i][j] = ru[i][grid_size_y-3];
    rvs[i][j] = rv[i][grid_size_y-3];
    energy_s[i][j] = energy[i][grid_size_y-3];
  }

  else
  {
    r[i][j] = rs[i][grid_size_y-3];
    ru[i][j] = rus[i][grid_size_y-3];
    rv[i][j] = rvs[i][grid_size_y-3];
    energy[i][j] = energy_s[i][grid_size_y-3];
  }
}

/*
  Corner point of rectangular object in path of flow. U and V should be initialized to 0 at startup.
  Pressure is taken to be an average of surrounding points.
*/
void MacCormack::BC_CORNER_POINT(size_t i, size_t j)
{
  if (temp[i][j] <= T_STACK) {
    temp[i][j] += HEAT_RATE;
  }
  p[i][j] = 0.25 * (p[i+1][j] + p[i-1][j] + p[i][j+1] + p[i][j-1]);

  if (predictor)
  {
    rs[i][j] = p[i][j] / (R*temp[i][j]);
    energy_s[i][j] = rs[i][j] * (cv * temp[i][j]);
  }

  else
  {
    r[i][j] = p[i][j] / (R*temp[i][j]);
    energy[i][j] = r[i][j] * (cv * temp[i][j]);
  }
}

void MacCormack::BC_RIGHT_PRESSURE_OUTLET(size_t i, size_t j)
{
  // extrapolating temperature and velocity
  temp[i][j] = 2.*temp[i-1][j] - temp[i-2][j];
  u[i][j] = 2.*u[i-1][j] - u[i-2][j];
  v[i][j] = 2.*v[i-1][j] - v[i-2][j];

  // temp[i][j] = (temp[i][j] + temp[i-1][j] + temp[i-2][j]) /3. ;
  // u[i][j] = (u[i][j] + u[i-1][j] + u[i-2][j]) /3.;
  // v[i][j] = (v[i][j] + v[i-1][j] + v[i-2][j]) /3.;

  // Using start-up initialized pressure to derive density update.

  if (predictor)
  {
    rs[i][j] = p[i][j] / (R*temp[i][j]);
    rus[i][j] = rs[i][j] * u[i][j];
    rvs[i][j] = rs[i][j] * v[i][j];
    energy_s[i][j] = rs[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  }

  else
  {
    r[i][j] = p[i][j] / (R*temp[i][j]);
    ru[i][j] = r[i][j] * u[i][j];
    rv[i][j] = r[i][j] * v[i][j];
    energy[i][j] = r[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  }

}

void MacCormack::BC_LEFT_INLET(size_t i, size_t j)
{
  /*
    extrapolating pressure for density update
  */
  p[i][j] = 2.*p[i+1][j] - p[i+2][j];

  /*
    Need to update momentums here in order to enforce a constant velocity constraint. Otherwise you will be 
    enforcing a constant momentum constraint. Not sure what the implications of a constant momentum constraint are. I
    guess that would be forcing a constant mass-flow at the inlet
  */

  if (predictor)
  {
    rs[i][j] = p[i][j] / (R*temp[i][j]);
    rus[i][j] = rs[i][j] * u[i][j];
    rvs[i][j] = 0.;
    energy_s[i][j] = rs[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  }

  else
  {
    r[i][j] = p[i][j] / (R*temp[i][j]);
    ru[i][j] = r[i][j] * u[i][j];
    rv[i][j] = 0.;
    energy[i][j] = r[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  }
}

void MacCormack::BC_RIGHT_OUTFLOW(size_t i, size_t j)
{
  if (predictor) 
  {
    rs[i][j] = 2.*r[i-1][j] - r[i-2][j];
    rus[i][j] = 2.*ru[i-1][j] - ru[i-2][j];
    rvs[i][j] = 2.*rv[i-1][j] - rv[i-2][j];
    energy_s[i][j] = 2.*energy[i-1][j] - energy[i-2][j];

    /*
      According to Anderson this is how the outflow should be done. Extrapolation of u, v, T, and P.
      This does not make sense to me since the stencil is solving for r, ru, rv, and E - then u, v, T, and P 
      are derived from those main quantities between predictor and corrector steps.
    */

    // u[i][j] = 2.*u[i-1][j] - u[i-2][j];
    // v[i][j] = 2.*v[i-1][j] - v[i-2][j];
    // temp[i][j] = 2.*temp[i-1][j] - temp[i-2][j];
    // p[i][j] = 2.*p[i-1][j] - p[i-2][j];

    // rs[i][j] = p[i][j] / (R*temp[i][j]);
    // rus[i][j] = rs[i][j] * u[i][j];
    // rvs[i][j] = rs[i][j] * v[i][j];
    // energy_s[i][j] = rs[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  } 
  
  else 
  {
    r[i][j] = 2.*rs[i-1][j] - rs[i-2][j];
    ru[i][j] = 2.*rus[i-1][j] - rus[i-2][j];
    rv[i][j] = 2.*rvs[i-1][j] - rvs[i-2][j];
    energy[i][j] = 2.*energy_s[i-1][j] - energy_s[i-2][j];

    // u[i][j] = 2.*u[i-1][j] - u[i-2][j];
    // v[i][j] = 2.*v[i-1][j] - v[i-2][j];
    // p[i][j] = 2.*p[i-1][j] - p[i-2][j];
    // temp[i][j] = 2.*temp[i-1][j] - temp[i-2][j];

    // r[i][j] = p[i][j] / (R*temp[i][j]);
    // ru[i][j] = r[i][j] * u[i][j];
    // rv[i][j] = r[i][j] * v[i][j];
    // energy[i][j] = r[i][j] * (cv * temp[i][j] + (u[i][j]*u[i][j] + v[i][j]*v[i][j]) / 2.);
  }
}

void MacCormack::BC_TOP_WALL(size_t i, size_t j)
{
  if (temp[i][j] <= T_STACK) {
    temp[i][j] += HEAT_RATE;
  }

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
  if (temp[i][j] <= T_STACK) {
    temp[i][j] += HEAT_RATE;
  }
  /* uncomment this for adiabatic wall condition */

  // temp[i][j] = temp[i][j+1];
  
  /*
    Anderson extrapolates pressure at the plate surface which can be used with temp to derive the density.
    This does result in a smoother density and pressure near the plate, although I doubt if it's physically valid.
    
    Using the usual one-sided difference for the density equation results in a jagged density pattern near the plate
    but this effect is diminished with higher grid resolution.
  */
  
  // p[i][j] = 2.*p[i][j+1] - p[i][j+2];

  if (predictor) {
    rs[i][j] = r[i][j] + dt/dy * 0.5 * (3.*rv[i][j] - 4.*rv[i][j+1] + rv[i][j+2]);
    // rs[i][j] = p[i][j] / (R*temp[i][j]);

    energy_s[i][j] = rs[i][j] * cv * temp[i][j];

  }

  else 
  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + dt/dy * 0.5 * (3.*rvs[i][j] - 4.*rvs[i][j+1] + rvs[i][j+2]));
    // r[i][j] = p[i][j] / (R*temp[i][j]);

    energy[i][j] = r[i][j] * cv * temp[i][j];
  }
}

void MacCormack::BC_RIGHT_WALL(size_t i, size_t j)
{ 
  if (temp[i][j] <= T_STACK) {
    temp[i][j] += HEAT_RATE;
  }
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
  if (temp[i][j] <= T_STACK) {
    temp[i][j] += HEAT_RATE;
  }
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

double MacCormack::calc_stencil(int component, size_t i, size_t j)
{
  if (component == 0) {
    return - dt/dx * (E0[i+leftx][j] - E0[i-rightx][j]) - dt/dy * (F0[i][j+upy] - F0[i][j-downy]);
  } else if (component == 1) {
    return - dt/dx * (E1[i+leftx][j] - E1[i-rightx][j]) - dt/dy * (F1[i][j+upy] - F1[i][j-downy]);
  } else if (component == 2) {
    return - dt/dx * (E2[i+leftx][j] - E2[i-rightx][j]) - dt/dy * (F2[i][j+upy] - F2[i][j-downy]);
  } else if (component == 3) {
    return - dt/dx * (E3[i+leftx][j] - E3[i-rightx][j]) - dt/dy * (F3[i][j+upy] - F3[i][j-downy]);
  }
}

void MacCormack::run_solver_step()
{
  if (TIMESTEP < 100) {
    HEAT_RATE = 0.;
  } else {
    HEAT_RATE = 0.00712;
  }

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
  update_E_and_F_Periodic();

  leftx = forward_diff_first;
  upy = forward_diff_first;
  rightx = !leftx;
  downy = !upy;
  
  // performing predictor step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == FREE_FLOW)
      { 
        rs[i][j] = r[i][j] + calc_stencil(0, i, j);
        rus[i][j] = ru[i][j] + calc_stencil(1, i, j);
        rvs[i][j] = rv[i][j] + calc_stencil(2, i, j);
        energy_s[i][j] = energy[i][j] + calc_stencil(3, i, j);
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

      if (region[i][j] != RIGHT_PRESSURE_OUTLET) {
        p[i][j] = rs[i][j] * R * temp[i][j];
      }
      
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * k_constant;
    }
  }

  predictor = false;
  update_tau_and_q();
  update_E_and_F();
  update_E_and_F_Periodic();

  // corrector runs opposite to predictor
  leftx = !forward_diff_first;
  upy = !forward_diff_first;
  rightx = !leftx;
  downy = !upy;

  // performing corrector step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {
      if (region[i][j] == FREE_FLOW)
      {       
        r[i][j] = 0.5 * (r[i][j] + rs[i][j] + calc_stencil(0, i, j));
        ru[i][j] = 0.5 * (ru[i][j] + rus[i][j] + calc_stencil(1, i, j));
        rv[i][j] = 0.5 * (rv[i][j] + rvs[i][j] + calc_stencil(2, i, j));
        energy[i][j] = 0.5 * (energy[i][j] + energy_s[i][j] + calc_stencil(3, i, j));
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

      if (region[i][j] != RIGHT_PRESSURE_OUTLET) {
        p[i][j] = r[i][j] * R * temp[i][j];
      }

      if (p[i][j] < 0.) {
        cout << "Negative pressure at " << i << " " << j << endl;
      }

      if (temp[i][j] < 0.) {
        cout << "Negative temperature at " << i << " " << j << endl;
      }

      if (r[i][j] < 0.) {
        cout << "Negative density at " << i << " " << j << endl;
      }
      
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * k_constant;
    }
  }

  cout << u[225][15] << endl;
  FILE * u_fp;
  u_fp = fopen("Data_Output/Probe.dat","a");
  fprintf(u_fp, "%.10lf ", p[450][15]);
  fclose(u_fp);

}
