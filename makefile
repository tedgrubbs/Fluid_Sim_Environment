LIBS=-fopenmp -L/usr/local/lib -lglfw3 -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp


Main: Main.o Simulation.o MacCormack.o glad.o
	ccache g++ -o Main Main.o Simulation.o MacCormack.o glad.o $(LIBS)

Simulation.o: Simulation.cpp Simulation.h
	ccache g++ -c  Simulation.cpp -o Simulation.o -fopenmp

MacCormack.o: MacCormack.cpp Simulation.h
	ccache g++ -c  MacCormack.cpp -o MacCormack.o -w -fopenmp

Main.o: Main.cpp Simulation.h
	ccache g++ -c  Main.cpp -o Main.o -fopenmp

glad.o: glad.c
	ccache g++ -c $(pkg-config --cflags glfw3) glad.c

clean:
	rm ./*.o

run:
	./Main
