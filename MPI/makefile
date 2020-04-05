OPT = -ffast-math -O3 -march=native
a.out: main.o functions.o
<<<<<<< HEAD
	mpicxx $(OPT) Zhordan_Main_MPI.o Zhordan_Functions_MPI.o  
=======
	mpicxx $(OPT) Zhordan_Main_MPI.cpp Zhordan_Functions_MPI.o  
>>>>>>> 50c29ce2add89106f2eabc0ab4a034750118eac8
main.o: Zhordan_Main_MPI.cpp Zhordan_Header_MPI.hpp
	mpicxx -c $(OPT) Zhordan_Main_MPI.cpp
functions.o: Zhordan_Functions_MPI.cpp Zhordan_Header_MPI.hpp
	mpicxx -c $(OPT) Zhordan_Functions_MPI.cpp
clean:
	rm -f Zhordan_Main_MPI.o Zhordan_Functions_MPI.o a.out
