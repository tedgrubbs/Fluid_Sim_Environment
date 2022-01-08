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
  dt = 0.05 * min_dim / max_boundary_speed;
  cout << "MacCormack timestep defined by stability criteria: " << dt << endl;

  rs = create2dArray<double>(grid_size_x, grid_size_y);
  us = create2dArray<double>(grid_size_x, grid_size_y);
  rus = create2dArray<double>(grid_size_x, grid_size_y);
  vs = create2dArray<double>(grid_size_x, grid_size_y);
  rvs = create2dArray<double>(grid_size_x, grid_size_y);
  ps = create2dArray<double>(grid_size_x, grid_size_y);
  temp_s = create2dArray<double>(grid_size_x, grid_size_y);
  energy_s = create2dArray<double>(grid_size_x, grid_size_y);
  int_energy_s = create2dArray<double>(grid_size_x, grid_size_y);

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

  a1 = dt / dx;
  a2 = dt / dy;
  a5 = 4. * dt  / (3.*dx*dx);
  a6 = dt / (dy*dy);
  a7 = dt / (dx*dx);
  a8 = 4.*dt / (3.*dy*dy);
  a9 = dt / (12. * dx*dy);
  a10 = 2.*(a5+a6);
  a11 = 2.*(a7+a8);

  b1 = 1./3.;
  b2 = 8. / (9.*dx*c*c);
  b3 = 1. / (18.*dy*c*c);
  b4 = 8.  / (9.*dy*c*c);
  b5 = 1. / (18.*dx*c*c);

  bool forward_diff_first = true;

}

/*
  General MacCormack formulation of compressible Navier-Stokes with simple alternation of difference direction. 
*/
void MacCormack::free_flow_predictor(size_t i, size_t j)
{

  int rightx,leftx,righty,lefty;
  if (forward_diff_first)
  {
    rightx = i+1;
    leftx = i;
    righty = j+1;
    lefty = j;
  }
  else
  {
    rightx = i;
    leftx = i-1;
    righty = j;
    lefty = j-1;
  }

  rs[i][j] = r[i][j] - a1 * (ru[rightx][j] - ru[leftx][j]) - a2 * (rv[i][righty] - rv[i][lefty]);

  rus[i][j] = ru[i][j] - a1 * (p[rightx][j] - p[leftx][j])
    - a1 * (r[rightx][j]*u[rightx][j]*u[rightx][j] - r[leftx][j]*u[leftx][j]*u[leftx][j])
    - a2 * (r[i][righty]*u[i][righty]*v[i][righty] - r[i][lefty]*u[i][lefty]*v[i][lefty])
    - mu[i][j] *  a10 * u[i][j]
    + mu[i][j] *a5 * (u[i+1][j] + u[i-1][j])
    +  mu[i][j] * a6 * (u[i][j+1] + u[i][j-1])
    + mu[i][j] *  a9 * (v[i+1][j+1] + v[i-1][j-1] - v[i+1][j-1] - v[i-1][j+1]);

  rvs[i][j] = rv[i][j] - a2 * (p[i][righty] - p[i][lefty])
    - a1 * (r[rightx][j]*u[rightx][j]*v[rightx][j] - r[leftx][j]*u[leftx][j]*v[leftx][j])
    - a2 * (r[i][righty]*v[i][righty]*v[i][righty] - r[i][lefty]*v[i][lefty]*v[i][lefty])
    -  mu[i][j] * a11 * v[i][j]
    +  mu[i][j] * a7 * (v[i+1][j] + v[i-1][j])
    +  mu[i][j] * a8 * (v[i][j+1] + v[i][j-1])
    + mu[i][j] *  a9 * (u[i+1][j+1] + u[i-1][j-1] - u[i+1][j-1] - u[i-1][j+1]);

  // energy derivatives
  double dxEPU = ((energy[rightx][j] + p[rightx][j]) * u[rightx][j] - (energy[leftx][j] + p[leftx][j]) * u[leftx][j]);
  double dyEPV = ((energy[i][righty] + p[i][righty]) * v[i][righty] - (energy[i][lefty] + p[i][lefty]) * v[i][lefty]);
  
  // u derivatives
  double dudx = (u[i+1][j] - u[i-1][j]) / (2.*dx);
  double du2dx2 = ((u[i+1][j] - 2.*u[i][j] + u[i-1][j]) / (dx*dx));
  double dudy = (u[i][j+1] - u[i][j-1]) / (2.*dy);
  double du2dy2 = ((u[i][j+1] - 2.*u[i][j] + u[i][j-1]) / (dy*dy));

  // v derivatives
  double dvdx = (v[i+1][j] - v[i-1][j]) / (2.*dx);
  double dv2dx2 = ((v[i+1][j] - 2.*v[i][j] + v[i-1][j]) / (dx*dx));
  double dvdy = (v[i][j+1] - v[i][j-1]) / (2.*dy);
  double dv2dy2 = ((v[i][j+1] - 2.*v[i][j] + v[i][j-1]) / (dy*dy));

  // mixed derivatives
  double du2dxdy = 1./(2.*dx) * ((u[i+1][j+1] - u[i+1][j-1])/ (2.*dy) - (u[i-1][j+1] - u[i-1][j-1]) / (2.*dy));
  double dv2dxdy = 1./(2.*dx) * ((v[i+1][j+1] - v[i+1][j-1])/ (2.*dy) - (v[i-1][j+1] - v[i-1][j-1]) / (2.*dy));

  // temperature derivatives
  double dT2dx2 = ((temp[i+1][j] - 2.*temp[i][j] + temp[i-1][j]) / (dx*dx));
  double dT2dy2 = ((temp[i][j+1] - 2.*temp[i][j] + temp[i][j-1]) / (dy*dy));

  energy_s[i][j] = energy[i][j] 
    - dt*(dxEPU - 2./3. * mu[i][j]  * (2.*(dudx*dudx + u[i][j]*du2dx2) - (dudx*dvdy + u[i][j]*dv2dxdy))
    - mu[i][j]  * ((dvdx*dudy + v[i][j]*du2dxdy) + (dvdx*dvdx + v[i][j]*dv2dx2)) - k[i][j] *dT2dx2)
    - dt*(dyEPV - 2./3. * mu[i][j]  * (2.*(dvdy*dvdy + v[i][j]*dv2dy2) - (dudx*dvdy + v[i][j]*du2dxdy))
    - mu[i][j]  * ((dudy*dudy + u[i][j]*du2dy2) + (dudy*dvdx + u[i][j]*dv2dxdy)) - k[i][j] *dT2dy2)
    ;

}

void MacCormack::free_flow_corrector(size_t i, size_t j)
{

  int rightx,leftx,righty,lefty;
  if (forward_diff_first)
  {
    rightx = i;
    leftx = i-1;
    righty = j;
    lefty = j-1;
  }
  else
  {
    rightx = i+1;
    leftx = i;
    righty = j+1;
    lefty = j;
  }

  r[i][j] = 0.5 * ((r[i][j] + rs[i][j])
    - a1 * (rus[rightx][j] - rus[leftx][j])
    - a2 * (rvs[i][righty] - rvs[i][lefty]));

  ru[i][j] = 0.5 * ((ru[i][j] + rus[i][j])
    - a1 * (ps[rightx][j] - ps[leftx][j])
    - a1 * (rs[rightx][j]*us[rightx][j]*us[rightx][j] - rs[leftx][j]*us[leftx][j]*us[leftx][j])
    - a2 * (rs[i][righty]*us[i][righty]*vs[i][righty] - rs[i][lefty]*us[i][lefty]*vs[i][lefty])
    - mu[i][j] *  a10 * us[i][j]
    + mu[i][j] * a5 * (us[i+1][j] + us[i-1][j])
    +  mu[i][j] * a6 * (us[i][j+1] + us[i][j-1])
    +  mu[i][j] * a9 * (vs[i+1][j+1] + vs[i-1][j-1] - vs[i+1][j-1] - vs[i-1][j+1]));

  rv[i][j] = 0.5 * ((rv[i][j] + rvs[i][j])
    - a2 * (ps[i][righty] - ps[i][lefty])
    - a1 * (rs[rightx][j]*us[rightx][j]*vs[rightx][j] - rs[leftx][j]*us[leftx][j]*vs[leftx][j])
    - a2 * (rs[i][righty]*vs[i][righty]*vs[i][righty] - rs[i][lefty]*vs[i][lefty]*vs[i][lefty])
    - mu[i][j] *  a11 * vs[i][j]
    +  mu[i][j] * a7 * (vs[i+1][j] + vs[i-1][j])
    +  mu[i][j] * a8 * (vs[i][j+1] + vs[i][j-1])
    + mu[i][j] *  a9 * (us[i+1][j+1] + us[i-1][j-1] - us[i+1][j-1] - us[i-1][j+1]));

  // energy derivatives
  double dxEPU = ((energy_s[rightx][j] + ps[rightx][j]) * us[rightx][j] - (energy_s[leftx][j] + ps[leftx][j]) * us[leftx][j]);
  double dyEPV =  ((energy_s[i][righty] + ps[i][righty]) * vs[i][righty] - (energy_s[i][lefty] + ps[i][lefty]) * vs[i][lefty]);
  
  // u derivatives
  double dudx = (us[i+1][j] - us[i-1][j]) / (2.*dx);
  double du2dx2 = ((us[i+1][j] - 2.*us[i][j] + us[i-1][j]) / (dx*dx));
  double dudy = (us[i][j+1] - us[i][j-1]) / (2.*dy);
  double du2dy2 = ((us[i][j+1] - 2.*us[i][j] + us[i][j-1]) / (dy*dy));

  // v derivatives
  double dvdx = (vs[i+1][j] - vs[i-1][j]) / (2.*dx);
  double dv2dx2 = ((vs[i+1][j] - 2.*vs[i][j] + vs[i-1][j]) / (dx*dx));
  double dvdy = (vs[i][j+1] - vs[i][j-1]) / (2.*dy);
  double dv2dy2 = ((vs[i][j+1] - 2.*vs[i][j] + vs[i][j-1]) / (dy*dy));

  // mixed derivatives
  double du2dxdy = 1./(2.*dx) * ((us[i+1][j+1] - us[i+1][j-1])/ (2.*dy) - (us[i-1][j+1] - us[i-1][j-1]) / (2.*dy));
  double dv2dxdy = 1./(2.*dx) * ((vs[i+1][j+1] - vs[i+1][j-1])/ (2.*dy) - (vs[i-1][j+1] - vs[i-1][j-1]) / (2.*dy));

  // temperature derivatives
  double dT2dx2 = ((temp[i+1][j] - 2.*temp[i][j] + temp[i-1][j]) / (dx*dx));
  double dT2dy2 = ((temp[i][j+1] - 2.*temp[i][j] + temp[i][j-1]) / (dy*dy));

  energy[i][j] = 0.5 * ((energy[i][j] + energy_s[i][j]) 
    - dt*(dxEPU - 2./3. * mu[i][j]  * (2.*(dudx*dudx + us[i][j]*du2dx2) - (dudx*dvdy + us[i][j]*dv2dxdy))
    - mu[i][j]  * ((dvdx*dudy + vs[i][j]*du2dxdy) + (dvdx*dvdx + vs[i][j]*dv2dx2)) - k[i][j] *dT2dx2)
    - dt*(dyEPV - 2./3. * mu[i][j]  * (2.*(dvdy*dvdy + vs[i][j]*dv2dy2) - (dudx*dvdy + vs[i][j]*du2dxdy))
    - mu[i][j]  * ((dudy*dudy + us[i][j]*du2dy2) + (dudy*dvdx + us[i][j]*dv2dxdy)) - k[i][j] *dT2dy2))
    ;
}

/*
  No-slip wall conditions for velocity, with density derived from continuity equation. Using 
  2nd-order accurate one-sided differences.
*/
void MacCormack::stationary_wall_predictor(size_t i, size_t j)
{

  rus[i][j] = 0.;
  rvs[i][j] = 0.;

  // stationary left wall
  if (region[i-1][j] == EXTERNAL)
  {
    rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  }

  // stationary right wall
  else if (region[i+1][j] == EXTERNAL)
  {
    rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
  }

  // stationary bottom wall
  else if (region[i][j-1] == EXTERNAL)
  {
    rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
  }

  // stationary top wall
  else if (region[i][j+1] == EXTERNAL)
  {
    rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
  }
  
  energy_s[i][j] = temp[i][j]*cv*rs[i][j];

}

void MacCormack::stationary_wall_corrector(size_t i, size_t j)
{

  ru[i][j] = 0.;
  rv[i][j] = 0.;

  // stationary left wall
  if (region[i-1][j] == EXTERNAL)
  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  }

  // stationary right wall
  else if (region[i+1][j] == EXTERNAL)
  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
  }

  // stationary bottom wall
  else if (region[i][j-1] == EXTERNAL)
  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
  }

  // stationary top wall
  else if (region[i][j+1] == EXTERNAL)
  {
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
  }

  energy[i][j] = temp[i][j]*cv*r[i][j];

}

// Moving wall/lid BC. Used when fluid should be moving parallel to surface, but not perpendicularly.
// Should take priority at corners when adjacent to stationary walls. Otherwise density will grow indefinitely at the corners
void MacCormack::moving_wall_predictor(size_t i, size_t j)
{

  // top moving lid
  if (region[i][j+1] == EXTERNAL)
  {
    // Check if this is a corner point
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID)
    {
      rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    }
    else
    {
      rs[i][j] = r[i][j] - 0.5*a1*boundary_v[i][j] * (r[i+1][j] - r[i-1][j]) + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    }

    // need this in case there is any y-gradient in velocity at the outlet
    if (region[i][j+1] == OUTLET || region[i][j-1] == OUTLET) rs[i][j] = 1.0;

    rus[i][j] = boundary_v[i][j] * rs[i][j];
    rvs[i][j] = 0.;
  }

  // bottom moving lid
  else if (region[i][j-1] == EXTERNAL)
  {
    // Check if this is a corner point
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID)
    {
      rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    }
    else
    {
      rs[i][j] = r[i][j] - 0.5*a1*boundary_v[i][j] * (r[i+1][j] - r[i-1][j]) - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    }

    if (region[i][j+1] == OUTLET || region[i][j-1] == OUTLET) rs[i][j] = 1.0;

    rus[i][j] = boundary_v[i][j] * rs[i][j];
    rvs[i][j] = 0.;
  }



}

void MacCormack::moving_wall_corrector(size_t i, size_t j)
{

  // top moving lid
  if (region[i][j+1] == EXTERNAL)
  {
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    }
    else
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*boundary_v[i][j] * (rs[i+1][j] - rs[i-1][j]) + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    }

    if (region[i][j+1] == OUTLET || region[i][j-1] == OUTLET) r[i][j] = 1.0;

    ru[i][j] = r[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
  }

  // bottom moving lid
  else if (region[i][j-1] == EXTERNAL)
  {
    if (region[i-1][j] != MOVING_LID || region[i+1][j] != MOVING_LID)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
    }
    else
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1*boundary_v[i][j] * (rs[i+1][j] - rs[i-1][j]) - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
    }

    if (region[i][j+1] == OUTLET || region[i][j-1] == OUTLET) r[i][j] = 1.0;

    ru[i][j] = r[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
  }


}

// velocity inlet condition with density varying according to the continuity equation.
// Should not take priority at outlet
void MacCormack::inlet_predictor(size_t i, size_t j)
{

  // inlet on left
  if (region[i-1][j] == EXTERNAL)
  {
    rus[i][j] = r[i][j] * boundary_v[i][j];
    rvs[i][j] = 0.;
    rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
  }

  // inlet on right
  else if (region[i+1][j] == EXTERNAL)
  {
    rus[i][j] = r[i][j] * boundary_v[i][j];
    rvs[i][j] = 0.;
    rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
  }

}

void MacCormack::inlet_corrector(size_t i, size_t j)
{

  // inlet on left
  if (region[i-1][j] == EXTERNAL)
  {
    ru[i][j] = rs[i][j] * boundary_v[i][j];
    rv[i][j] = 0.;
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
  }

  // inlet on right
  else if (region[i+1][j] == EXTERNAL)
  {
    rus[i][j] = r[i][j] * boundary_v[i][j];
    rvs[i][j] = 0.;
    rs[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
  }

}

/*
  Static boundary condition where all variables are specified, v is 0 but u can vary - hence the 'u' in the name
*/
void MacCormack::static_u_predictor(size_t i, size_t j) 
{
  rus[i][j] = r[i][j] * boundary_v[i][j];
  rvs[i][j] = 0.;
  rs[i][j] = r[i][j];
  energy_s[i][j] = rs[i][j]*(temp[i][j]*cv + .5*boundary_v[i][j]*boundary_v[i][j]) ;
}

void MacCormack::static_u_corrector(size_t i, size_t j) 
{
  ru[i][j] = rs[i][j] * boundary_v[i][j];
  rv[i][j] = 0.;
  r[i][j] = rs[i][j];
  energy[i][j] = r[i][j]*(temp_s[i][j]*cv + .5*boundary_v[i][j]*boundary_v[i][j]) ;
}

/*
  velocity outlet condition. Here only the U velocity is specified but density can vary according to the continuity equation.
  I observed more better looking results at higher reynolds numbers if I let v vary according to an average.
  the solution is still stable if you enforce a constant v, but a noticeable increase in velocity can be seen at the outlet due to
  the large gradient in y that occurs.

  should NOT take priority at corners
*/
void MacCormack::outlet_predictor(size_t i, size_t j)
{

  // outlet on right
  if (region[i+1][j] == EXTERNAL)
  {
    rus[i][j] = r[i][j] * boundary_v[i][j];
    rvs[i][j] = 1./4.*(rv[i-1][j] + rv[i][j+1] + rv[i][j-1] + rv[i][j]);
    rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]) - 0.5*a2*(rv[i][j+1] - rv[i][j-1]);
  }

}

void MacCormack::outlet_corrector(size_t i, size_t j)
{

  // outlet on right
  if (region[i+1][j] == EXTERNAL)
  {
    ru[i][j] = rs[i][j] * boundary_v[i][j];
    rv[i][j] = 1./4.*(rvs[i-1][j] + rvs[i][j+1] + rvs[i][j-1] + rvs[i][j]);
    r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]) - 0.5*a2*(rvs[i][j+1] - rvs[i][j-1]));
  }

}

/*
  extrapolated outlet conditions aka "outflow" boundaries. May cause simulation instability if used in compressible regimes 
  or when density or velocity varies dramatically out the outlet.
  Should take priority at the corner points
*/
void MacCormack::extrapolate_out_predictor(size_t i, size_t j)
{

  // outlet on right
  if (region[i+1][j] == EXTERNAL)
  {
    rus[i][j] = 2.*ru[i-1][j] - ru[i-2][j];
    rvs[i][j] = 2.*rv[i-1][j] - rv[i-2][j];
    rs[i][j] = 2.*r[i-1][j] - r[i-2][j];
    energy_s[i][j] = 2.*energy[i-1][j] - energy[i-2][j];

  }

}

void MacCormack::extrapolate_out_corrector(size_t i, size_t j)
{

  // outlet on right
  if (region[i+1][j] == EXTERNAL)
  {
    ru[i][j] = 2.*rus[i-1][j] - rus[i-2][j];
    rv[i][j] = 2.*rvs[i-1][j] - rvs[i-2][j];
    r[i][j] = 2.*rs[i-1][j] - rs[i-2][j];
    energy[i][j] = 2.*energy_s[i-1][j] - energy_s[i-2][j];
  }

}

/*
  No-slip conditions for stationary wall but using momentum equation to derive density. Only works if the equation of state is 
  p = rho * c^2.

  Defaults to continuity equation if at simulation corner.

  But if at corner of internal object - like the corner of a square - will use simple average of surrounding points to derive density.
*/
void MacCormack::stationary_wall_mom_predictor(size_t i, size_t j)
{

  rus[i][j] = 0.;
  rvs[i][j] = 0.;

  // if at a corner just default to normal stationary wall condition

  // left-facing wall (similar to the right wall)
  if (region[i+1][j] == EXTERNAL)
  {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL)
    {
      rs[i][j] = r[i][j] + 0.5*a1 * (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j]);
    }
    else
    {
      rs[i][j] = b1 * (4.*r[i-1][j] - r[i-2][j])
      + mu[i][j] *  b2 * (-5.*u[i-1][j] + 4.*u[i-2][j] - u[i-3][j])
      - mu[i][j] *  b3 * (
        -(v[i-2][j+1] - v[i-2][j-1])
        + 4.*(v[i-1][j+1] - v[i-1][j-1])
        -3.*(v[i][j+1] - v[i][j-1])
      );
    }
  }

  // right-facing wall (similar to the left wall)
  else if(region[i-1][j] == EXTERNAL)
  {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL)
    {
      rs[i][j] = r[i][j] - 0.5*a1 * (-ru[i+2][j] + 4.*ru[i+1][j] - 3.*ru[i][j]);
    }
    else
    {
      rs[i][j] = b1 * (4.*r[i+1][j] - r[i+2][j])
      -  mu[i][j] * b2 * (-5.*u[i+1][j] + 4.*u[i+2][j] - u[i+3][j])
      -  mu[i][j] * b3 * (
        -(v[i+2][j+1] - v[i+2][j-1])
        + 4.*(v[i+1][j+1] - v[i+1][j-1])
        -3.*(v[i][j+1] - v[i][j-1])
      );
    }
  }

  // upward-facing wall (similar to bottom wall)
  else if (region[i][j-1] == EXTERNAL)
  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL)
    {
      rs[i][j] = r[i][j] - 0.5*a2 * (-rv[i][j+2] + 4.*rv[i][j+1] - 3.*rv[i][j]);
    }
    else
    {
      rs[i][j] = b1 * (4.*r[i][j+1] - r[i][j+2])
      - mu[i][j] *  b4 * (-5.*v[i][j+1] + 4.*v[i][j+2] - v[i][j+3])
      - mu[i][j] *  b5 * (
        -(u[i+1][j+2] - u[i-1][j+2])
        + 4.*(u[i+1][j+1] - u[i-1][j+1])
        -3.*(u[i+1][j] - u[i-1][j])
      );
    }
  }

  // downward facing wall (similar to top wall)
  else if (region[i][j+1] == EXTERNAL)
  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL)
    {
      rs[i][j] = r[i][j] + 0.5*a2 * (-rv[i][j-2] + 4.*rv[i][j-1] - 3.*rv[i][j]);
    }
    else
    {
      rs[i][j] = b1 * (4.*r[i][j-1] - r[i][j-2])
      + mu[i][j] *  b4 * (-5.*v[i][j-1] + 4.*v[i][j-2] - v[i][j-3])
      - mu[i][j] *  b5 * (
        -(u[i+1][j-2] - u[i-1][j-2])
        + 4.*(u[i+1][j-1] - u[i-1][j-1])
        -3.*(u[i+1][j] - u[i-1][j])
      );
    }
  }
  else
  {
    rs[i][j] = .25 * (r[i+1][j] + r[i-1][j] + r[i][j+1] + r[i][j-1]);
  }

}

void MacCormack::stationary_wall_mom_corrector(size_t i, size_t j)
{

  ru[i][j] = 0.;
  rv[i][j] = 0.;

  // left-facing wall
  if (region[i+1][j] == EXTERNAL)
  {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a1 * (-rus[i-2][j] + 4.*rus[i-1][j] - 3.*rus[i][j]));
    }
    else
    {
      r[i][j] = b1 * (4.*rs[i-1][j] - rs[i-2][j])
      +  mu[i][j] * b2 * (-5.*us[i-1][j] + 4.*us[i-2][j] - us[i-3][j])
      - mu[i][j] *  b3 * (
        -(vs[i-2][j+1] - vs[i-2][j-1])
        + 4.*(vs[i-1][j+1] - vs[i-1][j-1])
        -3.*(vs[i][j+1] - vs[i][j-1])
      );
    }

  }

  // right-facing wall
  else if(region[i-1][j] == EXTERNAL)
  {
    if (region[i][j+1] == EXTERNAL || region[i][j-1] == EXTERNAL)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a1 * (-rus[i+2][j] + 4.*rus[i+1][j] - 3.*rus[i][j]));
    }
    else
    {
      r[i][j] = b1 * (4.*rs[i+1][j] - rs[i+2][j])
      - mu[i][j] *  b2 * (-5.*us[i+1][j] + 4.*us[i+2][j] - us[i+3][j])
      -  mu[i][j] * b3 * (
        -(vs[i+2][j+1] - vs[i+2][j-1])
        + 4.*(vs[i+1][j+1] - vs[i+1][j-1])
        -3.*(vs[i][j+1] - vs[i][j-1])
      );
    }

  }

  // upward-facing wall (similar to bottom wall)
  else if (region[i][j-1] == EXTERNAL)
  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] - 0.5*a2 * (-rvs[i][j+2] + 4.*rvs[i][j+1] - 3.*rvs[i][j]));
    }
    else
    {
      r[i][j] = b1 * (4.*rs[i][j+1] - rs[i][j+2])
      - mu[i][j] *  b4 * (-5.*vs[i][j+1] + 4.*vs[i][j+2] - vs[i][j+3])
      - mu[i][j] *  b5 * (
        -(us[i+1][j+2] - us[i-1][j+2])
        + 4.*(us[i+1][j+1] - us[i-1][j+1])
        -3.*(us[i+1][j] - us[i-1][j])
      );
    }
  }

  // downward facing wall (similar to top wall)
  else if (region[i][j+1] == EXTERNAL)
  {
    if (region[i-1][j] == EXTERNAL || region[i+1][j] == EXTERNAL)
    {
      r[i][j] = 0.5 * (r[i][j] + rs[i][j] + 0.5*a2 * (-rvs[i][j-2] + 4.*rvs[i][j-1] - 3.*rvs[i][j]));
    }
    else
    {
      r[i][j] = b1 * (4.*rs[i][j-1] - rs[i][j-2])
      + mu[i][j] *  b4 * (-5.*vs[i][j-1] + 4.*vs[i][j-2] - vs[i][j-3])
      - mu[i][j] *  b5 * (
        -(us[i+1][j-2] - us[i-1][j-2])
        + 4.*(us[i+1][j-1] - us[i-1][j-1])
        -3.*(us[i+1][j] - us[i-1][j])
      );
    }
  }
  else
  {
    r[i][j] = .25 * (rs[i+1][j] + rs[i-1][j] + rs[i][j+1] + rs[i][j-1]);
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

  // predictor step
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<(grid_size_x); ++i)
  {
    for (j=0; j<(grid_size_y); ++j)
    {

      if (region[i][j] == FREE_FLOW)
      {
        free_flow_predictor(i,j);
      }
      else if (region[i][j] == STATIONARY)
      {
        stationary_wall_predictor(i,j);
      }
      else if (region[i][j] == MOVING_LID)
      {
        moving_wall_predictor(i,j);
      }
      else if (region[i][j] == INLET)
      {
        inlet_predictor(i,j);
      }
      else if (region[i][j] == OUTLET)
      {
        outlet_predictor(i,j);
      }
      else if (region[i][j] == STATIONARY_MOM)
      {
        stationary_wall_mom_predictor(i,j);
      }
      else if (region[i][j] == STATIC_U)
      {
        static_u_predictor(i,j);
      }
      else if (region[i][j] == EX_OUTLET)
      {
        extrapolate_out_predictor(i,j);
      }


    }
  }

  // calculating starred velocities
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for ( i=0; i<grid_size_x; ++i)
  {
    for ( j=0; j<grid_size_y; ++j)
    {
      us[i][j] = rus[i][j] / rs[i][j];
      vs[i][j] = rvs[i][j] / rs[i][j];
      int_energy_s[i][j] = energy_s[i][j] / rs[i][j] - 0.5 * (us[i][j]*us[i][j] + vs[i][j]*vs[i][j]);
      temp[i][j] = int_energy_s[i][j] / cv;
      ps[i][j] = rs[i][j] * R * temp[i][j];
      // ps[i][j] = rs[i][j] *c*c;
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * cp / Pr;
      
    }
  }
  // cout << "Old pressure " << rs[32][32] *c*c << endl;
  // cout << "New k " << k[32][32] << endl;

  // Corrector step.
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for (i=0; i<grid_size_x; ++i)
  {
    for (j=0; j<grid_size_y; ++j)
    {
      if (region[i][j] == FREE_FLOW)
      {
        free_flow_corrector(i,j);
      }
      else if (region[i][j] == STATIONARY)
      {
        stationary_wall_corrector(i,j);
      }
      else if (region[i][j] == MOVING_LID)
      {
        moving_wall_corrector(i,j);
      }
      else if (region[i][j] == INLET)
      {
        inlet_corrector(i,j);
      }
      else if (region[i][j] == OUTLET)
      {
        outlet_corrector(i,j);
      }
      else if (region[i][j] == STATIONARY_MOM)
      {
        stationary_wall_mom_corrector(i,j);
      }
      else if (region[i][j] == STATIC_U)
      {
        static_u_corrector(i,j);
      }
      else if (region[i][j] == EX_OUTLET)
      {
        extrapolate_out_corrector(i,j);
      }

    }
  }

  // Calculating new velocties
  #pragma omp parallel for num_threads(MAX_THREADS) collapse(2) private(i,j)
  for ( i=0; i<grid_size_x; ++i)
  {
    for ( j=0; j<grid_size_y; ++j)
    {
      u[i][j] = ru[i][j] / r[i][j];
      v[i][j] = rv[i][j] / r[i][j];
      int_energy[i][j] = energy[i][j] / r[i][j]  - 0.5 * (u[i][j]*u[i][j] + v[i][j]*v[i][j]);
      temp[i][j] = int_energy[i][j] / cv;
      p[i][j] = r[i][j] * R * temp[i][j];
      // p[i][j] = r[i][j]* c*c;
      mu[i][j] = sutherland(temp[i][j]);
      k[i][j] = mu[i][j] * cp / Pr;
    }
  }

}
