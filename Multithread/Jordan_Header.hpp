#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <pthread.h>
#include <sys/time.h>
#define EPS 10e-17
using namespace std;
struct arg{
    int num;                    //номер потока
    int p;                      //количество потоков
    int n;                      //размер матрицы
    int m;                      //размер блока
    double *a;                  //ссылка на матрицу
    double *b;                  //ссылка на присоединенную матрицу
    char* file_name;            //имя файла
    double time;                //время работы
    int error;                  //ошибка
    pthread_barrier_t* barrier; //барреьр 
};
inline double matrixValue (int i, int j);
double discrepansy(double *oldMatrix, double *matrix, int n,int m, double*row);
double blockNorm(double* matrix, int rows, int columns);
double getIJ(double* matrix, int i, int j, int n, int m);
double matrixNorm(double* matrix, int n, int m);
double discrepansy2(double* matrix, double* inv_matrix,  int n, int m);
void getBlock(const double* matrix, double* block, int g, int h, int n, int m);
void putBlock(double* matrix, const double* block, int g, int h, int n, int m);
void createIdMatrix(double* matrix, int n, int m);
void createMatrix(double* matrix, int n, int m);
void printMatrix(const double* matrix, int n, int m, int num_print);
void matrixMulti(double* A, int a1, int a2, double* B, int b2, double* Res);
void printBlock(double *block, int n, int m);
void blockMulti(double* A, double* B, double* Res, int n, int m);
void createBlock(double* block, int n, int m);
void blockSubsruction(double* A, double* B, int n, int m);
void swapColumn(double* block, int i, int j, int n);
void swapRows(double* block, int i, int j, int n);
void multiRow(double* block, double c, int i, int j, int n, int m);
void matrixId(double* block, int n);
void swapMatrixColumn(double* matrix, int i, int j, int n, int m, int whole, int remainder);
void swapRowsInArray(int* temp2, int i, int j);
void swapMatrixRows(double* matrix, int i, int j, int n, int m, int whole, int remainder);
void matrixThreadsZeros(arg* a);
void matrixThreadsInverse(arg* a, double matrix_norm, double* block, double* help1, int* temp, int* temp2, double* inv_block);
void findMainElement(arg* a, int i, double norm, int* index, double* min_norm, double* temp, double* inv_block, double* block);
void matrixThreadsNorm(arg* a, double* norm);
void* thread_func(void* args);
int blockInverse(double* block, int* temp2,  double* id, int n, double norm);
int readMatrixFromFile(double* matrix, int n, int m, string file_name);
int matrixInverse(double* matrix, double* res, int n, int m, int* temp, int* temp2, double* block, double*  help1, double* inv_block);
int matrixThreadsData(arg* a);