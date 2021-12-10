#include "Simulation.h"


std::chrono::high_resolution_clock::time_point begin,end;
std::chrono::microseconds duration;

void update(int value) {



}

int main(int argc, char const *argv[]) {

  begin = std::chrono::high_resolution_clock::now();
  MacCormack sim;
  sim.run();

  return 0;
}
