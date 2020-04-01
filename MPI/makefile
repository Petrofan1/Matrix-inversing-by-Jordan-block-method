OPT = -ffast-math -O3 -march=native
a.out: main.o functions.o
	mpicxx $(OPT) Zhordan_Main_MPI.cpp Zhordan_Functions_MPI.o  
main.o: Zhordan_Main_MPI.cpp Zhordan_Header_MPI.hpp
	mpicxx -c $(OPT) Zhordan_Main_MPI.cpp
functions.o: Zhordan_Functions_MPI.cpp Zhordan_Header_MPI.hpp
	mpicxx -c $(OPT) Zhordan_Functions_MPI.cpp
clean:
	rm -f main.o functions.o a.out
