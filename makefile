Main: Main.o
	ccache g++ -o Main Main.o -fopenmp

Main.o: Main.cpp
	ccache g++ -c  Main.cpp -o Main.o -fopenmp

clean:
	rm Main.o

run:
	./Main
