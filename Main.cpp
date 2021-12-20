#include "Simulation.h"

int main(int argc, char const *argv[]) 
{
  std::chrono::high_resolution_clock::time_point begin,end;
  std::chrono::microseconds duration;
  begin = std::chrono::high_resolution_clock::now();

  MacCormack sim;
  sim.run();

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
  std::cout << "Runtime: " << duration.count()/1000. << " milliseconds" << std::endl;

  return 0;
}
