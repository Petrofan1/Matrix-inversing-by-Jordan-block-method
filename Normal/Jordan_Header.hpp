#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#define EPS 10e-18
using namespace std;

void getBlock(const double* matrix, double* block, int g, int h, int n, int m);
void putBlock(double* matrix, const double* block, int g, int h, int n, int m);
double blockNorm(double* matrix, int rows, int columns);
inline double matrixValue (int i, int j);
int readMatrixFromFile(double* matrix, int n, int m, string file_name);
void createIdMatrix(double* matrix, int n, int m);
void createMatrix(double* matrix, int n, int m);
void printMatrix(const double* matrix, int n, int m, int num_print);
void matrixMulti(double* A, int a1, int a2, double* B, int b2, double* Res);
double getIJ(double* matrix, int i, int j, int n, int m);
void printBlock(double *block, int n, int m);
void blockMulti(double* A, double* B, double* Res, int n, int m);
void createBlock(double *block, int n, int m);
void blockSubsruction(double* A, double* B, int n, int m);
void swapColumn(double* block, int i, int j, int n);
void swapRows(double* block, int i, int j, int n);
void multiRow(double* block, double c, int i, int j, int n, int m);
void matrixId(double* block, int n);
void swapMatrixColumn(double* matrix, int i, int j, int n, int m, int whole, int remainder);
void swapRowsInArray(int* temp2, int i, int j);
void swapMatrixRows(double* matrix, int i, int j, int n, int m, int whole, int remainder);
int blockInverse(double* block, int* temp2,  double* id, int n, double norm);
double matrixNorm(double* matrix, int n, int m);
int matrixInverse(double* matrix, double* res, int n, int m, int* temp, int* temp2, double* block, double*  help1, double* inv_block);
double discrepansy2(double* matrix, double* inv_matrix,  int n, int m);