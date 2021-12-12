LIBS=-fopenmp -L/usr/local/lib -lglfw3 -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp
COMPILEFLAGS=-c -fopenmp 

Main: Main.o Simulation.o MacCormack.o glad.o
	ccache g++ -o Main Main.o Simulation.o MacCormack.o glad.o $(LIBS)

Simulation.o: Simulation.cpp Simulation.h
	ccache g++ $(COMPILEFLAGS) Simulation.cpp -o Simulation.o 

MacCormack.o: MacCormack.cpp Simulation.h
	ccache g++ $(COMPILEFLAGS)  MacCormack.cpp -o MacCormack.o

Main.o: Main.cpp Simulation.h
	ccache g++ $(COMPILEFLAGS) Main.cpp -o Main.o

glad.o: glad.c
	ccache g++ -c $(pkg-config --cflags glfw3) glad.c

clean:
	rm ./*.o

run:
	./Main
