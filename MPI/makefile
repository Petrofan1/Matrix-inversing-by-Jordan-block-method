OPT = -ffast-math -O3 -march=native
a.out: main.o functions.o
	mpicxx $(OPT) ZhordanMainMPI.cpp ZhordanFunctionsMPI.o  
main.o: ZhordanMainMPI.cpp ZhordanHeaderMPI.hpp
	mpicxx -c $(OPT) ZhordanMainMPI.cpp
functions.o: ZhordanFunctionsMPI.cpp ZhordanHeaderMPI.hpp
	mpicxx -c $(OPT) ZhordanFunctionsMPI.cpp
clean:
	rm -f main.o functions.o a.out
