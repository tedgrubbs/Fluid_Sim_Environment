LIBS=-fopenmp -L/usr/local/lib -lglfw3 -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp


Main: Main.o Utilities.o Solver.o glad.o
	ccache g++ -o Main Main.o Utilities.o Solver.o glad.o $(LIBS)

Utilities.o: Utilities.cpp Utilities.h
	ccache g++ -c  Utilities.cpp -o Utilities.o -fopenmp

Solver.o: Solver.cpp Solver.h Utilities.h
	ccache g++ -c  Solver.cpp -o Solver.o -w -fopenmp

Main.o: Main.cpp Utilities.h Solver.h
	ccache g++ -c  Main.cpp -o Main.o -fopenmp

glad.o: glad.c
	ccache g++ -c $(pkg-config --cflags glfw3) glad.c

clean:
	rm Main.o Solver.o Utilities.o glad.o

run:
	./Main
