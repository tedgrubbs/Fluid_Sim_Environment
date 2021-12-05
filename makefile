LIBS=-lglut -lGL -fopenmp


Main: Main.o Utilities.o Solver.o
	ccache g++ -o Main Main.o Utilities.o Solver.o $(LIBS)

Utilities.o: Utilities.cpp Utilities.h
	ccache g++ -c  Utilities.cpp -o Utilities.o -fopenmp

Solver.o: Solver.cpp Solver.h Utilities.h
	ccache g++ -c  Solver.cpp -o Solver.o -w -fopenmp

Main.o: Main.cpp Utilities.h Solver.h
	ccache g++ -c  Main.cpp -o Main.o -fopenmp

clean:
	rm Main.o Solver.o Utilities.o

run:
	./Main
