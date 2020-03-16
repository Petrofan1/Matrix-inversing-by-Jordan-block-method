OPT =  -ffast-math -O3 -march=native -pthread
a.out: main.o functions.o
	g++ $(OPT) Jordan_Main.cpp Jordan_Functions_Parallel.o  
main.o: Jordan_Main.cpp Jordan_Header.hpp 
	g++ -c $(OPT) Jordan_Main.cpp
functions.o: Jordan_Functions_Parallel.cpp Jordan_Header.hpp
	g++ -c $(OPT) Jordan_Functions_Parallel.cpp
clean:
	rm -f main.o functions.o a.out
