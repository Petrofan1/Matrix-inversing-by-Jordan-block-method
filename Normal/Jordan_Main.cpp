
#include"Jordan_Header.hpp"
using namespace std;
int main(int argc, char* argv[])
	/*  argv[1] = matrix_size -- размер матрицы;
	    argv[2] = block_size -- размер блока;
	    Значение argc должно быть не меньше трех: a.out + n + m; */
{
	if (argc < 3)
	{
		cout<<"Недостаточно параметров для выплнения программы."<<endl;
		return 1;
	}
	if (atoi(argv[1]) < atoi(argv[2]))
	{
		cout<<"Введенный размер блока больше размера матрицы."<<endl;
		return 1;
	}

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
        double *matrix1, *matrix2, *matrix3, *block, *inv_block, *help1, *temp3;
	int *temp2, *temp;
	double t;
	matrix3 = new double[n*n];
	matrix1 = new double[n*n];
	matrix2 = new double[n*n];
	block = new double[m*m];
	inv_block = new double[m*m];
	help1 = new double[m*m];
	temp = new int[m*m];
	temp2 = new int[n/m];
        temp3 = new double[m];
	for(int i = 0; i < n/m; i++)
	{
		temp2[i]=i;
	}
	cout<<"Создаю матрицы..."<<endl;
	if(argc == 4)
	{
		if(readMatrixFromFile(matrix1, n, m, argv[3]) == -1)
		{
			delete[] matrix1;
			delete[] matrix2;
			delete[] matrix3;
			delete[] block;
			delete[] help1;
			delete[] inv_block;
			delete[] temp;
			delete[] temp2;
                    delete[] temp3;
			return -1;
		};
	}
	else
	{
		createMatrix(matrix1, n, m);
	}
	memcpy(matrix3, matrix1, sizeof(double)*n*n);
	cout<<"Обращаю матрицу..."<<endl;
	t=clock();
        if(n<=8)
        {
            printMatrix(matrix1, n,m,n);
        }
        else
        {
            printMatrix(matrix1, n,m,8);
        }
	createIdMatrix(matrix2, n, m);
        //printMatrix(matrix1, n, m, n);
	if(matrixInverse(matrix1, matrix2, n, m, temp, temp2, block, help1, inv_block)==-1)
	{
		delete[] matrix1;
		delete[] matrix2;
		delete[] matrix3;
		delete[] block;
				delete[] help1;
				delete[] inv_block;
				delete[] temp;
				delete[] temp2;
            delete[] temp3;
		return -1;
        }
        //printMatrix(matrix2, n, m, n);
        cout<<"\nTime = "<<(clock()-t)/CLOCKS_PER_SEC<<"\n"<<endl;
        cout<<"Residual = "<<discrepansy2(matrix3,n,m, matrix2, temp3)<<endl;
	delete[] matrix1;
	delete[] matrix2;
	delete[] matrix3;
	delete[] block;
        delete[] help1;
        delete[] inv_block;
        delete[] temp;
        delete[] temp2;
        delete[] temp3;
	return 0;
}
