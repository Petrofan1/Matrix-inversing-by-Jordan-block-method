OPT = -ffast-math -O3 -march=native -pg
a.out: main.o functions.o
	g++ $(OPT) Jordan_Main.o Jordan_Functions.o  
main.o: Jordan_Main.cpp Jordan_Header.hpp
	g++ -c $(OPT) Jordan_Main.cpp
functions.o: Jordan_Functions.cpp Jordan_Header.hpp
	g++ -c $(OPT) Jordan_Functions.cpp
clean:
	rm -f Jordan_Main.o Jordan_Functions.o a.out
