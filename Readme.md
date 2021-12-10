# Fluid_Sim_V2

The code here is meant to serve as a simulation environment for fluid mechanics problems. The code is designed in a way that
allows for new CFD methods to be added without much copying and pasting of the same old code.

A base Simulation class is used to store all variables that should be present in all solvers. Solver classes inherit the Simulation class and must define the virtual run_solver_step() method in order to be valid. An example here is given with a simple MacCormack solver. One can see that in the MacCormack.cpp file the new class simply defines extra variables that are specific to the MacCormack scheme and then defines the run_solver_step() method. All other all simulation functions- graphics, data recording, stop conditions, etc - is handled by the base Simulation class.

### Note:
On my systems it is definitely slower to implement the solver function as a class method. A faster option would be to define a solver function that operates only on a globally defined struct or object that contains all relevant simulation variables. This was done in my first attempt of implementing the MacCormack scheme and it is indeed faster than the currently implemented object-oriented (OO) code. However, I believe that the OO style implemented here allows the user to implement new schemes that are easier to debug and understand. If more speed is required, one can always return to global variables.
