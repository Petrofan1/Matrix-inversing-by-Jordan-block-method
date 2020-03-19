#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iomanip>
#include <stdio.h>
//#include "mpich/mpi.h"
#include "mpi.h"

using namespace std;

#define SETW_CONSTANT 12	//количество печатаемых символов для элементов матрицы
#define EPS 1e-15               //эпсилон
#define PRINT_CONSTANT 10	//размера печатаемой части матрицы
#define TYPE_OF_MATRIX 1   	//матричная формула: 	0 -- единичная
                                //                      1 -- |i - j|
                                //                      2 -- Гильберт
struct buffer{
    double value;
    int index;
};

double getTime();
inline double matrixValue (int i, int j, int mode);
double blockNorm(double* block, int n, int m);
void printBlock(double* block, int n, int m);
void printMatrix(const double* matrix, int n, int m);
void blockSubsruction(double* A, double* B, int n, int m);
void blockSum(double* A, double* B, int n, int m);
void matrixMulti(double* A, int a1, int a2, double* B, int b2, double* Res);
void matrixId(double* block, int n);
void swapColumn(double* block, int i, int j, int n);
void swapRows(double* block, int i, int j, int n);
void swapRowsInArray(int* array, int i, int j);
int blockInverse(double* block, int* temp2,  double* id, int n, double norm);
int readMatrixFromFile(double* matrix, double* row, int n, int m, int p, int my_rank, char* file_name);
void createMatrix(double* matrix, int n, int m, int p, int my_rank, int type);
void MPINorm(double* matrix, double* row, int n, int m, int my_rank, int p, double* norm);
void printMatrixMPI(double* matrix, int my_rank, int p, int n, int m, double* row);
int matrixInverseMPI(double* matrix, double* id, double* row_a, double* row_b, double norm, int n, int m, int p, int my_rank);
void findMinElementMPI(int i, int n, int m, int p, int my_rank, double norm, int* main_index, double* main_value, double* row, int* swap_in_block, double* inv_block, double* block);
void swapMatrixColumnMPI(int n, int m, int p, int my_rank, int i, int j, double* matrix);
void findIndexOfMainElement(double* invec, double* inoutvec, int* len, MPI_Datatype* dtype);
double residualMPI();


