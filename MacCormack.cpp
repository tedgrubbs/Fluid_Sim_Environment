#include "Simulation.h"

size_t i;
size_t j;

MacCormack::MacCormack() : Simulation() {};

void MacCormack::run_solver_step() {
}
//
// void run_sim_timestep() {
//
//   // predictor step for interior points
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2) private(i,j)
//   for (i=1; i<(sim.grid_size_x-1); ++i) {
//     for (j=1; j<(sim.grid_size_y-1); ++j) {
//
//       sim.rs[i][j] = sim.r[i][j] - sim.a1 * (sim.ru[i+1][j] - sim.ru[i][j]) - sim.a2 * (sim.rv[i][j+1] - sim.rv[i][j]);
//
//       sim.rus[i][j] = sim.ru[i][j] - sim.a3 * (sim.r[i+1][j] - sim.r[i][j])
//         - sim.a1 * (sim.r[i+1][j]*sim.u[i+1][j]*sim.u[i+1][j] - sim.r[i][j]*sim.u[i][j]*sim.u[i][j])
//         - sim.a2 * (sim.r[i][j+1]*sim.u[i][j+1]*sim.v[i][j+1] - sim.r[i][j]*sim.u[i][j]*sim.v[i][j])
//         - sim.a10 * sim.u[i][j]
//         + sim.a5 * (sim.u[i+1][j] + sim.u[i-1][j])
//         + sim.a6 * (sim.u[i][j+1] + sim.u[i][j-1])
//         + sim.a9 * (sim.v[i+1][j+1] + sim.v[i-1][j-1] - sim.v[i+1][j-1] - sim.v[i-1][j+1]);
//
//       sim.rvs[i][j] = sim.rv[i][j] - sim.a4 * (sim.r[i][j+1] - sim.r[i][j])
//         - sim.a1 * (sim.r[i+1][j]*sim.u[i+1][j]*sim.v[i+1][j] - sim.r[i][j]*sim.u[i][j]*sim.v[i][j])
//         - sim.a2 * (sim.r[i][j+1]*sim.v[i][j+1]*sim.v[i][j+1] - sim.r[i][j]*sim.v[i][j]*sim.v[i][j])
//         - sim.a11 * sim.v[i][j]
//         + sim.a7 * (sim.v[i+1][j] + sim.v[i-1][j])
//         + sim.a8 * (sim.v[i][j+1] + sim.v[i][j-1])
//         + sim.a9 * (sim.u[i+1][j+1] + sim.u[i-1][j-1] - sim.u[i+1][j-1] - sim.u[i-1][j+1]);
//
//     }
//   }
//
//   // left and right wall boundary conditions
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(1) private(i,j)
//   for (j=0; j<sim.grid_size_y; ++j) {
//
//     i = 0;
//     sim.rs[i][j] = sim.r[i][j] - 0.5*sim.a1 * (-sim.ru[i+2][j] + 4.*sim.ru[i+1][j] - 3.*sim.ru[i][j]);
//     // The nature of the upwind difference at the lid edges causes density at the corners to increase indefinitely.
//     // I think the density at these corners should not differ too greatly from the neighboring gridpoints so
//     // doing an average instead of the original difference makes more physical sense.
//     if (j == (sim.grid_size_y-1)) {
//       sim.rs[i][j] = 1./3. * (sim.r[i][j] + sim.r[i+1][j] + sim.r[i][j-1]);
//     }
//
//     sim.rus[i][j] = 0.;
//     sim.rvs[i][j] = 0.;
//
//     i = (sim.grid_size_x-1);
//     sim.rs[i][j] = sim.r[i][j] + 0.5*sim.a1 * (-sim.ru[i-2][j] + 4.*sim.ru[i-1][j] - 3.*sim.ru[i][j]);
//     if (j == (sim.grid_size_y-1)) {
//       sim.rs[i][j] = 1./3. * (sim.r[i][j] + sim.r[i-1][j] + sim.r[i][j-1]);
//     }
//
//     sim.rus[i][j] = 0.;
//     sim.rvs[i][j] = 0.;
//
//   }
//
//   // top and bottom boundary conditions
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(1) private(i,j)
//   for (i=1; i<(sim.grid_size_x-1); ++i) {
//     j = 0;
//     sim.rs[i][j] = sim.r[i][j] - 0.5*sim.a2 * (-sim.rv[i][j+2] + 4.*sim.rv[i][j+1] - 3.*sim.rv[i][j]);
//     sim.rus[i][j] = 0.;
//     sim.rvs[i][j] = 0.;
//
//     j = (sim.grid_size_x-1);
//     sim.rs[i][j] = sim.r[i][j] - 0.5*sim.a1*sim.u_lid * (sim.r[i+1][j] - sim.r[i-1][j]) + 0.5*sim.a2 * (-sim.rv[i][j-2] + 4.*sim.rv[i][j-1] - 3.*sim.rv[i][j]);
//
//     sim.rus[i][j] = sim.u_lid * sim.rs[i][j];
//     sim.rvs[i][j] = 0.;
//   }
//
//   // corrector step for interior points
//
//   // calculating starred velocities
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2) private(i,j)
//   for ( i=0; i<sim.grid_size_x; ++i) {
//     for ( j=0; j<sim.grid_size_y; ++j) {
//       sim.us[i][j] = sim.rus[i][j] / sim.rs[i][j];
//       sim.vs[i][j] = sim.rvs[i][j] / sim.rs[i][j];
//     }
//   }
//
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2) private(i,j)
//   for (i=1; i<(sim.grid_size_x-1); ++i) {
//     for (j=1; j<(sim.grid_size_y-1); ++j) {
//
//       sim.r[i][j] = 0.5 * ((sim.r[i][j] + sim.rs[i][j])
//         - sim.a1 * (sim.rus[i][j] - sim.rus[i-1][j])
//         - sim.a2 * (sim.rvs[i][j] - sim.rvs[i][j-1]));
//
//       sim.ru[i][j] = 0.5 * ((sim.ru[i][j] + sim.rus[i][j])
//         - sim.a3 * (sim.rs[i][j] - sim.rs[i-1][j])
//         - sim.a1 * (sim.rs[i][j]*sim.us[i][j]*sim.us[i][j] - sim.rs[i-1][j]*sim.us[i-1][j]*sim.us[i-1][j])
//         - sim.a2 * (sim.rs[i][j]*sim.us[i][j]*sim.vs[i][j] - sim.rs[i][j-1]*sim.us[i][j-1]*sim.vs[i][j-1])
//         - sim.a10 * sim.us[i][j]
//         + sim.a5 * (sim.us[i+1][j] + sim.us[i-1][j])
//         + sim.a6 * (sim.us[i][j+1] + sim.us[i][j-1])
//         + sim.a9 * (sim.vs[i+1][j+1] + sim.vs[i-1][j-1] - sim.vs[i+1][j-1] - sim.vs[i-1][j+1]));
//
//       sim.rv[i][j] = 0.5 * ((sim.rv[i][j] + sim.rvs[i][j])
//         - sim.a4 * (sim.rs[i][j] - sim.rs[i][j-1])
//         - sim.a1 * (sim.rs[i][j]*sim.us[i][j]*sim.vs[i][j] - sim.rs[i-1][j]*sim.us[i-1][j]*sim.vs[i-1][j])
//         - sim.a2 * (sim.rs[i][j]*sim.vs[i][j]*sim.vs[i][j] - sim.rs[i][j-1]*sim.vs[i][j-1]*sim.vs[i][j-1])
//         - sim.a11 * sim.vs[i][j]
//         + sim.a7 * (sim.vs[i+1][j] + sim.vs[i-1][j])
//         + sim.a8 * (sim.vs[i][j+1] + sim.vs[i][j-1])
//         + sim.a9 * (sim.us[i+1][j+1] + sim.us[i-1][j-1] - sim.us[i+1][j-1] - sim.us[i-1][j+1]));
//
//     }
//   }
//
//   // left and right wall boundary conditions
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(1) private(i,j)
//   for (j=0; j<sim.grid_size_y; ++j) {
//     i = 0;
//     sim.r[i][j] = 0.5 * (sim.r[i][j] + sim.rs[i][j] - 0.5*sim.a1 * (-sim.rus[i+2][j] + 4.*sim.rus[i+1][j] - 3.*sim.rus[i][j]));
//     // need to add trick for corners, otherwise you get infinite densities here
//     if (j == (sim.grid_size_y-1)) {
//       sim.r[i][j] = 1./3. * (sim.rs[i][j] + sim.rs[i+1][j] + sim.rs[i][j-1]);
//     }
//
//     sim.ru[i][j] = 0.;
//     sim.rv[i][j] = 0.;
//
//     i = (sim.grid_size_x-1);
//     sim.r[i][j] = 0.5 * (sim.r[i][j] + sim.rs[i][j] + 0.5*sim.a1 * (-sim.rus[i-2][j] + 4.*sim.rus[i-1][j] - 3.*sim.rus[i][j]));
//     if (j == (sim.grid_size_y-1)) {
//       sim.r[i][j] = 1./3. * (sim.rs[i][j] + sim.rs[i-1][j] + sim.rs[i][j-1]);
//     }
//
//     sim.ru[i][j] = 0.;
//     sim.rv[i][j] = 0.;
//   }
//
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(1) private(i,j)
//   for (i=1; i<(sim.grid_size_x-1); ++i) {
//     j = 0;
//     sim.r[i][j] = 0.5 * (sim.r[i][j] + sim.rs[i][j] - 0.5*sim.a2 * (-sim.rvs[i][j+2] + 4.*sim.rvs[i][j+1] - 3.*sim.rvs[i][j]));
//     sim.ru[i][j] = 0.;
//     sim.rv[i][j] = 0.;
//
//     j = (sim.grid_size_y-1);
//     sim.r[i][j] = 0.5 * (sim.r[i][j] + sim.rs[i][j] - 0.5*sim.a1*sim.u_lid * (sim.rs[i+1][j] - sim.rs[i-1][j]) + 0.5*sim.a2 * (-sim.rvs[i][j-2] + 4.*sim.rvs[i][j-1] - 3.*sim.rvs[i][j]));
//
//     sim.ru[i][j] = sim.r[i][j] * sim.u_lid;
//     sim.rv[i][j] = 0.;
//   }
//
//
//   #pragma omp parallel for num_threads(info_struct.MAX_THREADS) collapse(2) private(i,j)
//   for ( i=0; i<sim.grid_size_x; ++i) {
//     for ( j=0; j<sim.grid_size_y; ++j) {
//       sim.u[i][j] = sim.ru[i][j] / sim.r[i][j];
//       sim.v[i][j] = sim.rv[i][j] / sim.r[i][j];
//     }
//   }
//
// }
