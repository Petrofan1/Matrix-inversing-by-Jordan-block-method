#include "Zhordan_Header_MPI.hpp"
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

double getTime()
{
	struct timeval t;
	gettimeofday(&t, 0);
	return (double)t.tv_sec+(double)t.tv_usec/1000000.;
}
inline double matrixValue (int i, int j, int type)
{
	if(type == 0)
	{
		if(i == j) return 1;
		else return 0;
	}
	if(type == 1)
	{
		return fabs(i-j);
	}
	if(type == 2)
	{
		return 1/(double(i + j + 1));
	}
	return 0;
}
double blockNorm(double* block, int n, int m)
{
	double max = 0;
	double current;
	for(int i = 0; i < n; i++)
	{
		current = 0;
		for(int j = 0; j < m; j++)
		{
			current += fabs(block[i * n + j]);
		}
		if(current > max)
		{
			max = current;
		}
	}
	return max;
}
void MPINorm(double* matrix, double* row, int n, int m, int my_rank, int p, double* norm)
{
	int current, j;
	int whole = n/m, remainder = n % m, tag = 0;
	memset(row, 0, __SIZEOF_DOUBLE__*n);
	int counter = 0;
	for(j = my_rank ; j < whole; j+=p)
	{
		for(int k = 0; k < m; k++)
		{
			current = 0;
			while(current < whole)
			{
				for(int i = 0; i < m; i++)
				{
					row[j] += fabs(matrix[i + current*m*m + k*m + counter*n*m]);
				}
				current++;
			}
			if(remainder != 0)
			{
				for(int i = 0; i < remainder; i++)
				{
					row[j] += fabs(matrix[i + whole*m*m + k*remainder + counter*n*m]);
				}
			}
		}
		counter++;
	}
	if(remainder != 0 && (whole - my_rank) % p == 0)
	{
		for(int k = 0; k < remainder; k++)
		{
			current = 0;
			while(current < whole)
			{
				for(int i = 0; i < m; i++)
				{
					row[j] += fabs(matrix[i + current*m*remainder + k*m + counter*n*m]);
				}
				current++;
			}
			for(int i = 0; i < remainder; i++)
			{
				row[j] += fabs(matrix[i + whole*m*remainder + k*remainder + counter*n*m]);
			}
		}
	}
	if(my_rank != 0) MPI_Send (row, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	else
	{
		double* buf = new double[n];
		MPI_Status status;
		double max;
		for(int j = 1; j < p; j++)
		{
			MPI_Recv(buf, n, MPI_DOUBLE, MPI_ANY_SOURCE , tag, MPI_COMM_WORLD, &status);
			for(int i = 0; i < n; i++) row[i] += buf[i];
		}
		max = row[0];
		for(int i = 1; i < n; i++) max = (max > row[i]) ? max: row[i];
		*norm = max;
		delete []  buf;
	}
        MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}
void printBlock(double* block, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			cout<<setw(SETW_CONSTANT)<<block[j + i*m]<<" ";
		}
		cout<<endl;
	}
	cout<<endl;
}
void blockSubsruction(double* A, double* B, int n, int m)
{
	for(int i = 0; i < n*m; i++)
	{
		A[i] -= B[i];
	}
}
void blockSum(double* A, double* B, int n, int m)
{
	for(int i = 0; i < n*m; i++)
	{
		A[i] += B[i];
	}
}
void blockMinusId(double* block, int m)
{
	for(int i = 0; i < m; i++)
	{
		block[i*m + i] -= 1.0;
	}
}
void matrixMulti(double* A, int a1, int a2, double* B, int b2, double* Res)
{
        double *pB=B, *pA=A, *pRes=Res, sum[9];
        int remainder_i = a1 % 3;
        int remainder_j = b2 % 3;
        double a0, a9, a3, b0, b1, b3;
        for(int i = 0; i < a1 - remainder_i; i+=3)
        {
                pRes = Res + i*b2;
                for(int j = 0; j < b2 - remainder_j; j+=3)
                {
                        pA = A + i*a2;
                        pB = B  + j;
                        memset(sum, 0, sizeof(double)*9);
                        for(int k = 0; k < a2; k++)
                        {
                                a0 = pA[0];
                                a9 = pA[a2];
                                a3 = pA[2*a2];
                                b0 = pB[0];
                                b1 = pB[1];
                                b3 = pB[2];
                                sum[0] += a0*b0;
                                sum[1] += a0*b1;
                                sum[2] += a0*b3;
                                sum[3] += a9*b0;
                                sum[4] += a9*b1;
                                sum[5] += a9*b3;
                                sum[6] += a3*b0;
                                sum[7] += a3*b1;
                                sum[8] += a3*b3;
                                pA++;
                                pB+=b2;
                        }
                        pRes[0] = sum[0];
                        pRes[1] = sum[1];
                        pRes[2] = sum[2];
                        pRes[b2] = sum[3];
                        pRes[b2 + 1] = sum[4];
                        pRes[b2 + 2] = sum[5];
                        pRes[2*b2] = sum[6];
                        pRes[2*b2 + 1] = sum[7];
                        pRes[2*b2 + 2] = sum[8];
                        pRes += 3;
                }
        }
        for(int i = a1 - remainder_i; i < a1; i++)
        {
                for(int j = b2 - remainder_j; j < b2; j++)
                {
                        sum[0] = 0.;
                        pRes = Res + i*b2 + j;
                        for (int k = 0; k < a2; k++)
                        {
                                pA = A + i*a2 + k;
                                pB = B + k*b2 + j;
                                sum[0] += pA[0]*pB[0];
                        }
                        pRes[0] = sum[0];
                }
        }
        for(int i = 0; i < a1 - remainder_i; i++)
        {
                for(int j = b2 - remainder_j; j < b2; j++)
                {
                        sum[0] = 0.;
                        pRes = Res + i*b2 + j;
                        for (int k = 0; k < a2; k++)
                        {
                                pA=A + i*a2 + k;
                                pB=B + k*b2 + j;
                                sum[0] += pA[0]*pB[0];
                        }
                        pRes[0] = sum[0];
                }
        }
        for(int i = a1 - remainder_i; i < a1; i++)
        {
                for(int j = 0; j < b2 - remainder_j; j++)
                {
                        sum[0]=0.;
                        pRes = Res + i*b2 + j;
                        for(int k = 0; k < a2; k++)
                        {
                                pA = A + i*a2 + k;
                                pB = B + k*b2 + j;
                                sum[0] += pA[0]*pB[0];
                        }
                        pRes[0] = sum[0];
                }
        }
}
void matrixId(double* block, int n) 
{
	for(int s = 0; s < n; s++)
	{
		for(int p = 0; p < n; p++)
		{
			block[s*n + p] = (s == p) ? 1 : 0;
		}
	}
}
void swapColumn(double* block, int i, int j, int n)
{
	double temp;
	for(int k = 0; k < n; k++)
	{
		temp = block[k*n + j];
		block[k*n + j] = block[k*n + i];
		block[k*n + i] = temp;
	}
}
void swapRows(double* block, int i, int j, int n)
{
	double temp;
	for(int k = 0; k < n; k++)
	{
		temp = block[i*n + k];
		block[i*n + k] = block [j*n + k];
		block[j*n + k] = temp;
	}
}
void swapRowsInArray(int* array, int i, int j)
{
	int temp;
	temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}
void swapMatrixColumnMPI(int n, int m, int p, int my_rank, int i, int j, double* matrix, double* block)
{
	int remainder = n % m;
	int whole = n/m;
	int counter = 0;
	int nm = n*m;
	int mm = m*m;
	int mr = m*remainder;
	size_t size_m = m*m*(__SIZEOF_DOUBLE__);
	size_t size_r = m*remainder*(__SIZEOF_DOUBLE__);
	for (int k = my_rank;  k < whole; k += p)
	{
		memcpy(block, matrix + (counter*nm + j*mm), size_m);
		memcpy(matrix + (counter*nm + j*mm), matrix + (counter*nm + i*mm), size_m);
		memcpy(matrix + (counter*nm + i*mm), block,  size_m);
		counter++;
	}
	if((whole - my_rank) % p == 0 && remainder != 0)
	{
		memcpy(block, matrix + (counter*nm + j*mr), size_r);
		memcpy(matrix + (counter*nm + j*mr), matrix + (counter*nm + i*mr), size_r);
		memcpy(matrix + (counter*nm + i*mr), block,  size_r);
	}
}
void swapRowsInMatrixMPI(int n, int m, int p, int my_rank, int i, int j, double* row_a, double* row_b, double* matrix)
{
	int tag = 0;
	MPI_Status status;
        int temp = 0;
	if(i % p == my_rank && j % p == my_rank)
	{
		memcpy(row_a, matrix + (i/p)*n*m, __SIZEOF_DOUBLE__*n*m);
		memcpy(matrix + (i/p)*n*m, matrix + (j/p)*n*m, __SIZEOF_DOUBLE__*n*m);
		memcpy(matrix + (j/p)*n*m, row_a, __SIZEOF_DOUBLE__*n*m);
	}
	else 
	{
		if(i % p == my_rank)
		{
			memcpy(row_a, matrix + (i/p)*n*m, __SIZEOF_DOUBLE__*n*m);
			MPI_Sendrecv(row_a, n*m, MPI_DOUBLE, j % p, tag, row_b, n*m, MPI_DOUBLE, j % p, tag, MPI_COMM_WORLD, &status);
			memcpy(matrix + (i/p)*n*m, row_b, __SIZEOF_DOUBLE__*n*m);
		}
		if(j % p == my_rank)
		{
			memcpy(row_a, matrix + (j/p)*n*m, __SIZEOF_DOUBLE__*n*m);
			MPI_Sendrecv(row_a, n*m, MPI_DOUBLE, i % p, tag, row_b, n*m, MPI_DOUBLE, i % p, tag, MPI_COMM_WORLD, &status);
			memcpy(matrix + (j/p)*n*m, row_b, __SIZEOF_DOUBLE__*n*m);
		}
	}
        MPI_Bcast(&temp, 1, MPI_INT, i % p, MPI_COMM_WORLD);
        MPI_Bcast(&temp, 1, MPI_INT, j % p, MPI_COMM_WORLD);
}
int blockInverse(double* block, int* temp2,  double* id, int n, double norm)
{
	int max;
	int s;
	double MAX;
	for(int i = 0; i < n; i++)
	{
		temp2[i] = i;
	}
	matrixId(id, n);
	for(int j = 0; j < n; j++)
	{
		max = j;
		MAX = 0.0;
		for(int i = j; i < n; i++)
		{
			if(MAX < fabs(block[j*n + i]))
			{
				max = i;
				MAX = fabs(block[j*n + i]);
			}
		}
		if(MAX < EPS*norm)
		{
			return -1;
		}
		swapColumn(block, j, max, n);
		s = temp2[j];
		temp2[j] = temp2[max];
		temp2[max] = s;
		for(int k = j+1; k < n; k++)
		{
			block[j*n + k] /=(block[j*n + j]);
		}
		for(int k = 0; k < n; k++)
		{
			id[j*n + k] /=(block[j*n + j]);
		}
		for(int l = 0; l < j; l++)
		{
			for(int k = j + 1; k < n; k++)
			{
				block[l*n + k] -= (block[j*n + k])*block[l*n + j];
			}
			for(int k = 0; k < n; k++)
			{
				id[l*n + k] -= (id[j*n + k])*block[l*n + j];
			}
		}
		for(int l = j + 1; l < n; l++)
		{
			for(int k = j + 1; k < n; k++)
			{
				block[l*n + k] -= (block[j*n + k])*block[l*n + j];
			}
			for(int k = 0; k < n; k++)
			{
				id[l*n + k] -= (id[j*n + k])*block[l*n + j];
			}
		}

	}
	for(int l = 0; l < n; )
	{
		if(temp2[l] != l)
		{
			swapRows(id, temp2[l], l, n);
			swapRowsInArray(temp2, temp2[l], l);
		}
		else l++;
	}
	return 0;
}
int readMatrixFromFileMPI(double* matrix, double* row, int n, int m, int p, int my_rank, char* file_name)
{
	FILE *file;
	MPI_Status status;
	double current, temp;
	int counter = 0, remainder, whole, old_row = 0, new_row = 0, pos = 0, amount = 1, block_row = 0, tag = 0, l;
	int error = 0;
	remainder = n % m;
	whole = n/m;
	if(my_rank == 0)
	{
		file = fopen(file_name, "r");       		//Предполагается, что матрица в файле ровно того размера,
		if (file == NULL)                           //что был указан при запуске программы. Этот факт дополнительно не проверяется.
		{
			cout<<"\nCannot open file\n"<<endl;
			error = -1;
			MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
			return -1;
		}
		MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
		while(fscanf(file, "%lf", &temp) == 1) counter++;
		if(counter != n*n)
		{
				cout<<"\nНе хватает данных\n"<<endl;
				error = -1;
				MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
				return -1;
		}
		MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
		counter = 0;
		rewind(file);
		while(fscanf(file, "%lf", &temp) == 1)
		{
			counter++;
			if(old_row % m == 0 && old_row > 0)
			{
				new_row++;
				pos = 0;
				amount = 1;
				old_row = 0;
				row[pos] = temp;
				amount++; 
				pos++;
				continue;
			}
			row[pos] = temp;
			if(counter % (n*m) == 0)
			{
				if(block_row % p == 0)
				{
					memcpy(matrix + n*m*block_row/p , row , sizeof(double)*m*n);
				}
				else
				{
					MPI_Send(row, n*m, MPI_DOUBLE, block_row % p, tag, MPI_COMM_WORLD);
				}
				block_row++;
			}
			if(remainder != 0 && counter == n*n)
			{
				if(block_row % p == 0) memcpy(matrix + n*m*block_row/p , row , sizeof(double)*remainder*n);
				else
				{
					MPI_Send(row, remainder*n, MPI_DOUBLE, block_row % p, tag, MPI_COMM_WORLD);
				}
			}
			if(amount == n)
			{
				old_row++;
				amount=1;
				pos = m * old_row;
				continue;
			}
			if(amount % m == 0)
			{
				if((amount + remainder) != n || remainder == 0)
				{
					pos = amount*m + old_row*m;
				}
				else
				{
					pos = amount*m + remainder*old_row;
				}
				if(new_row*m + remainder == n)
				{
					pos = amount*remainder + old_row*m;
				}
				if(new_row*m + remainder == n && (amount + remainder) == n)
				{
					pos = amount*remainder + old_row*remainder;
				}
				amount++;
				continue;
			}
			amount++;
			pos++;
		}
		fclose(file);
	}
	else
	{
		current = my_rank;
		l = 0;
                MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if(error == -1) return -1;
                MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if(error == -1) return -1;
		while (current < whole)
		{
			MPI_Recv(matrix + l*m*n, m*n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
			l++;
			current += p;
		}
		if(remainder != 0)
		{
			MPI_Recv(matrix + l*m*n, remainder*n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
		} 
	}
	return 0;
}
void createMatrixFormula(double* matrix, int n, int m, int p, int my_rank, int type)
{
	int current = my_rank;
	int whole = n/m;
	int remainder = n % m;
	int l = 0;
	while(current < whole)
	{
		for(int k = 0; k < whole; k++)
		{
			for(int i = 0; i < m; i++)
			{
				for(int j = 0; j < m; j++)
				{
					matrix[j + i*m + k*m*m + l*n*m] = matrixValue(i + current*m, j + k*m, type);
				}

			}
		}
		if(remainder != 0)
		{
			for(int i = 0; i < m; i++)
			{
				for(int j = 0; j < remainder; j++)
				{
					matrix[j + i*remainder + whole*m*m + l*n*m] = matrixValue(i + current*m, j + whole*m, type);
				}
			}
		}
		l++;
		current += p;
	}
	if(remainder != 0 && whole % p  == my_rank)
	{
		for(int k = 0; k < whole; k++)
		{
			for(int i = 0; i < remainder; i++)
			{
				for(int j = 0; j < m; j++)
				{
					matrix[j + i*m + k*m*remainder + l*n*m] = matrixValue(i + current*m, j + k*m, type);
				}
			}
		}
		for(int i = 0; i < remainder; i++)
		{
			for(int j = 0; j < remainder; j++)
			{
				matrix[j + i*remainder + whole*m*remainder + l*n*m] = matrixValue(i + current*m, j + whole*m, type);
			}

		}

	}
}
int createMatrix(double* matrix, double* row, int n, int m, int p, int my_rank, int file_or_formula, char* file_name)
{
    if(file_or_formula == 4)
    {
        if(readMatrixFromFileMPI(matrix, row, n, m, p, my_rank, file_name) == -1)
        {
            return -1;
        }
    }
    else
    {
        createMatrixFormula(matrix, n, m, p, my_rank, TYPE_OF_MATRIX);
    }
    return 0;
}
void printMatrixMPI(double* matrix, int my_rank, int p, int n, int m, double* row)
{
	int whole = n/m;
	int remainder = n % m;
	int tag = 0;
	int current = my_rank;
	int l = 0;
	int LOCAL_PRINT_CONSTANT = (((int)PRINT_CONSTANT > n) ? n: (int)PRINT_CONSTANT);
	int print_v = 0;
	int print_h = 0;
	MPI_Status status;
	if(my_rank == 0)
	{
		for(int i = 0; i < whole; i++)
		{
			if(i % p == 0)
			{
				memcpy(row, matrix + (i/p)*n*m, sizeof(double)*n*m);
			}
			else
			{
				MPI_Recv(row, m*n, MPI_DOUBLE, i % p, tag, MPI_COMM_WORLD, &status);
			}
			for(int l = 0; l < m; l++)
			{
				print_h = 0;
				for(int k = 0; k < whole; k++)
				{
					for(int j = 0; j < m; j++)
					{
						cout<<setw(SETW_CONSTANT)<<row[j + k*m*m + l*m];
						print_h++;
						if(print_h == LOCAL_PRINT_CONSTANT) break;
					}
					if(print_h == LOCAL_PRINT_CONSTANT) break;
				}
				if(remainder != 0 && print_h != LOCAL_PRINT_CONSTANT)
				{
					for(int j = 0; j < remainder; j++)
					{
						cout<<setw(SETW_CONSTANT)<<row[j + whole*m*m + l*remainder];
						print_h++;
						if(print_h == LOCAL_PRINT_CONSTANT) break;
					}
				}
				printf("\n");
				print_v++;
				if(print_v == LOCAL_PRINT_CONSTANT) break;
			}
			if(print_v == LOCAL_PRINT_CONSTANT) break;
		}
		if(remainder != 0 && print_v != LOCAL_PRINT_CONSTANT)
		{
			if(whole % p == 0)
			{
				memcpy(row, matrix + (whole/p)*n*m, sizeof(double)*n*remainder);
			}
			else
			{
				MPI_Recv(row, remainder*n, MPI_DOUBLE, whole % p , tag, MPI_COMM_WORLD, &status);
			}
			for(int l = 0; l < remainder; l++)
			{
				print_h = 0;
				for(int k = 0; k < whole; k++)
				{
					for(int j = 0; j < m; j++)
					{
						cout<<setw(SETW_CONSTANT)<<row[j + k*m*remainder + l*m];
						print_h++;
						if(print_h == LOCAL_PRINT_CONSTANT) break;
					}
					if(print_h == LOCAL_PRINT_CONSTANT) break;
				}
				if(remainder != 0 && print_h != LOCAL_PRINT_CONSTANT)
				{
					for(int j = 0; j < remainder; j++)
					{
						cout<<setw(SETW_CONSTANT)<<row[j + whole*remainder*m + l*remainder];
						print_h++;
						if(print_h == LOCAL_PRINT_CONSTANT) break;
					}
				}
				printf("\n");
				print_v++;
				if(print_v == LOCAL_PRINT_CONSTANT) break;
			}
		}
		cout<<endl;
	}
	else
	{
		while (current < whole && current*m < LOCAL_PRINT_CONSTANT)
		{
			MPI_Send(matrix + l*m*n, n*m, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			current += p;
			l++;
		}
		if(remainder != 0 && whole % p == my_rank && current*m < LOCAL_PRINT_CONSTANT)
		{
			MPI_Send(matrix + l*m*n, n*remainder, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		}
	}
}
void findIndexOfMainElement(double* invec, double* inoutvec, int* len, MPI_Datatype *dtype)
{
	cout<<"Find index"<<endl;
	if(invec[0] > 0)
	{
		if(inoutvec[0] > invec[0] || inoutvec[0] < 0)
		{
			inoutvec[0] = invec[0];
			inoutvec[1] = invec[1];
		}
	}
	len = len;
	dtype = dtype;
}
void findMinElementMPI(int i, int n, int m, int p, int my_rank, double norm, int* main_index, double* main_value, double* row, int* swap_in_block, double* inv_block, double* block)
{
	int whole = n/m;
	int current_min_index = -1;
	double current_min_value = DBL_MAX;
	double start = my_rank;
	double block_norm = 0;
	while(start < i)
	{
		start += p;
	}
	for(int j = start; j < whole; j += p)
	{
		memcpy(block, row + j*m*m, m*m*__SIZEOF_DOUBLE__);
		if(blockInverse(block, swap_in_block, inv_block, m, norm) == -1)
		{
			continue;
		}
		block_norm = blockNorm(inv_block, m, m);
		if(block_norm < current_min_value)
		{
			current_min_index = j;
			current_min_value = block_norm;
		}
	}
	*main_index = current_min_index;
	*main_value = current_min_value;
}
int matrixInverseMPI(double* matrix, double* id, double* row_a, double* row_b, double norm, int n, int m, int p, int my_rank)
{
	int remainder = n % m;
	int whole = n/m;
	int nm = n*m;
	int mm = m*m;
	int mr = m*remainder;
	int rm = mr;
	int rr = remainder*remainder;
	int nr = n*remainder;
	int main_index, s, counter;
	int* swap_in_matrix;
	int* swap_in_block;
	double main_value;
	double* block;
	double* inv_block;
	size_t size_m = m*m*(__SIZEOF_DOUBLE__);
	size_t size_r = m*remainder*(__SIZEOF_DOUBLE__);
	swap_in_matrix = new int[whole];
	swap_in_block = new int[m];
	block = new double[m*m];
	inv_block = new double[m*m];
	buffer buffer_in;
	buffer buffer_out;
	buffer_in.value = DBL_MAX;
	buffer_in.index = -1;
	buffer_out.value = DBL_MAX;
	buffer_out.index = -1;
	for(int i = 0; i < whole; i++)
	{
		swap_in_matrix[i]=i;
	}
	if(remainder != 0)
	{
		for(int i = 0; i < whole; i++)
		{
			if(i % p == my_rank)
			{
				memcpy(row_a, matrix + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
				memcpy(row_b, id + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
			}
			MPI_Bcast(row_a, nm, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			MPI_Bcast(row_b, nm, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			findMinElementMPI(i, n, m, p, my_rank, norm, &main_index, &main_value, row_a, swap_in_block, inv_block, block);
			buffer_in.value = main_value;
			buffer_in.index = main_index;
			buffer_out.value = -1;
			buffer_out.index = -1;
			MPI_Allreduce(&buffer_in, &buffer_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
			main_index = buffer_out.index;
			if(main_index == -1)
			{
				delete [] swap_in_block;
				delete [] swap_in_matrix;
				delete [] block;
				delete [] inv_block;
				return -1;
			}
			memcpy(block, row_a + main_index*mm, size_m);
			if(main_index != i)
			{
				swapMatrixColumnMPI(n, m, p, my_rank, i, main_index, matrix, block);
				s = swap_in_matrix[i];
				swap_in_matrix[i] = swap_in_matrix[main_index];
				swap_in_matrix[main_index] = s;
				memcpy(row_a + main_index*mm, row_a + i*mm, size_m);
				memcpy(row_a + i*mm, block, size_m);
			}
			blockInverse(block, swap_in_block, inv_block, m, norm);
			if(i % p == my_rank)
			{
				for(int k = i + 1; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(inv_block, m, m, matrix + mm*k + (i/p)*nm, m, block);
					memcpy(matrix + mm*k + (i/p)*nm, block, size_m);
				}
				memset(block, 0, size_m);
				matrixMulti(inv_block, m, m, matrix + whole*mm + (i/p)*nm, remainder, block);
				memcpy(matrix + mm*whole + (i/p)*nm, block, size_r);
				for(int k = 0; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(inv_block, m, m, id + mm*k + (i/p)*nm, m, block);
					memcpy(id + mm*k + (i/p)*nm, block, size_m);
				}
				memset(block, 0, size_m);
				matrixMulti(inv_block, m, m, id + (i/p)*nm + whole*mm, remainder, block);
				memcpy(id + (i/p)*nm + mm*whole, block, size_r);
				memcpy(row_a, matrix + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
				memcpy(row_b, id + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
			}
			else
			{
				for(int k = i + 1; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(inv_block, m, m, row_a + mm*k, m, block);
					memcpy(row_a + mm*k, block, size_m);
				}
				memset(block, 0, size_m);
				matrixMulti(inv_block, m, m, row_a + whole*mm, remainder, block);
				memcpy(row_a + mm*whole, block, size_r);
				for(int k = 0; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(inv_block, m, m, row_b + mm*k, m, block);
					memcpy(row_b + mm*k, block, size_m);
				}
				memset(block, 0, size_m);
				matrixMulti(inv_block, m, m, row_b + whole*mm, remainder, block);
				memcpy(row_b + mm*whole, block, size_r);
			}
			counter = 0;
			for(int l = my_rank; l < whole; l += p)
			{
				if(l != i)
				{
					for(int k = i + 1; k < whole; k++)
					{
						memset(block, 0, size_m);
						matrixMulti(matrix + (counter*nm + i*mm), m, m, row_a + k*mm, m, block);
						blockSubsruction(matrix + (counter*nm + k*mm), block, m, m);
					}
					memset(block, 0, size_m);
					matrixMulti(matrix + (counter*nm + i*mm), m, m, row_a + whole*mm, remainder, block);
					blockSubsruction(matrix + (counter*nm + whole*mm), block, m, remainder);
					for(int k = 0; k < whole; k++)
					{
						memset(block, 0, size_m);
						matrixMulti(matrix + (counter*nm + i*mm), m, m, row_b + k*mm, m, block);
						blockSubsruction(id + (counter*nm + k*mm), block, m, m);
					}
					memset(block, 0, size_m);
					matrixMulti(matrix + (counter*nm + i*mm), m, m, row_b + whole*mm, remainder, block);
					blockSubsruction(id + (counter*nm + whole*mm), block, m, remainder);
				}
				counter++;
			}
			if((whole - my_rank) % p == 0)
			{
				for(int k = i + 1; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(matrix + (counter*nm + i*rm), remainder, m, row_a + k*mm,  m, block);
					blockSubsruction(matrix + (counter*nm + k*rm), block, remainder, m);
				}
				memset(block, 0, size_m);
				matrixMulti(matrix + (counter*nm + i*rm), remainder, m, row_a + whole*mm, remainder, block);
				blockSubsruction(matrix + (counter*nm + whole*rm), block, remainder, remainder);
				for(int k = 0; k < whole; k++)
				{
					memset(block, 0, size_m);
					matrixMulti(matrix + (counter*nm + i*rm), remainder, m, row_b + k*mm, m, block);
					blockSubsruction(id + (counter*nm + k*rm), block, remainder, m);
				}
				memset(block, 0, size_m);
				matrixMulti(matrix + (counter*nm + i*rm), remainder, m, row_b + whole*mm, remainder, block);
				blockSubsruction(id + (counter*nm + whole*rm), block, remainder, remainder);
			}
		}
		if(whole % p == my_rank)
		{
			memcpy(row_a, matrix + (whole/p)*nm, __SIZEOF_DOUBLE__*nm);
			memcpy(row_b, id + (whole/p)*nm, __SIZEOF_DOUBLE__*nm);
		}
		MPI_Bcast(row_a, nr, MPI_DOUBLE, whole % p, MPI_COMM_WORLD);
		MPI_Bcast(row_b, nr, MPI_DOUBLE, whole % p, MPI_COMM_WORLD);
		memcpy(block, row_a + whole*rm, __SIZEOF_DOUBLE__*rr);
		if(blockInverse(block, swap_in_block, inv_block, remainder, norm) == -1)
		{
			delete [] swap_in_block;
			delete [] swap_in_matrix;
			delete [] block;
			delete [] inv_block;
			return -1;
		}
		if(whole % p == my_rank)
		{
			for(int k = 0; k < whole; k++)
			{
				memset(block, 0, size_m);
				matrixMulti(inv_block, remainder, remainder, id + rm*k + (whole/p)*nm, m, block);
				memcpy(id + (whole/p)*nm + k*rm, block, size_r);
			}
			memset(block, 0, size_m);
			matrixMulti(inv_block, remainder, remainder, id + (whole/p)*nm + whole*rm,  remainder, block);
			memcpy(id + (whole/p)*nm + whole*rm, block, __SIZEOF_DOUBLE__*rr);
			memcpy(row_a, matrix + (whole/p)*nm, __SIZEOF_DOUBLE__*nm);
			memcpy(row_b, id + (whole/p)*nm, __SIZEOF_DOUBLE__*nm);
		}
		else
		{
			for(int k = 0; k < whole; k++)
			{
				memset(block, 0, size_m);
				matrixMulti(inv_block, remainder, remainder, row_b + rm*k, m, block);
				memcpy(row_b + k*rm, block, size_r);
			}
			memset(block, 0, size_m);
			matrixMulti(inv_block, remainder, remainder, row_b + whole*rm,  remainder, block);
			memcpy(row_b + whole*rm, block, __SIZEOF_DOUBLE__*rr);
		}
		counter = 0;
		for(int l = my_rank; l < whole; l += p)
		{
			for(int k = 0; k < whole; k++)
			{
				memset(block, 0, size_m);
				matrixMulti(matrix + (counter*nm + whole*mm), m, remainder, row_b + k*rm,  m, block);
				blockSubsruction(id + (counter*nm + k*mm), block, m, m);
			}
			memset(block, 0, size_m);
			matrixMulti(matrix + (counter*nm + whole*mm), m, remainder, row_b + whole*rm, remainder, block);
			blockSubsruction(id + (counter*nm + whole*mm), block, m, remainder);
			counter++;
		}
	}
	else
	{
		for(int i = 0; i < whole; i++)
		{
			if(i % p == my_rank)
			{
				memcpy(row_a, matrix + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
				memcpy(row_b, id + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
			}
			MPI_Bcast(row_a, nm, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			MPI_Bcast(row_b, nm, MPI_DOUBLE, i % p, MPI_COMM_WORLD);
			findMinElementMPI(i, n, m, p, my_rank, norm, &main_index, &main_value, row_a, swap_in_block, inv_block, block);
			buffer_in.value = main_value;
			buffer_in.index = main_index;
			buffer_out.value = -1;
			buffer_out.index = -1;
			MPI_Allreduce(&buffer_in, &buffer_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
			main_index = buffer_out.index;
			if(main_index == -1)
			{
				delete [] swap_in_block;
				delete [] swap_in_matrix;
				delete [] block;
				delete [] inv_block;
				return -1;
			}
			memset(block, 0, size_m);
			memcpy(block, row_a + main_index*mm, size_m);
			if(main_index != i)
			{
				swapMatrixColumnMPI(n, m, p, my_rank, i, main_index, matrix, inv_block);
				s = swap_in_matrix[i];
				swap_in_matrix[i] = swap_in_matrix[main_index];
				swap_in_matrix[main_index] = s;
				memcpy(row_a + main_index*mm, row_a + i*mm, size_m);
				memcpy(row_a + i*mm, block, size_m);
			}
			blockInverse(block, swap_in_block, inv_block, m, norm);
			if(i % p == my_rank)
			{
				for(int k = i + 1; k < whole; k++)
				{
					matrixMulti(inv_block, m, m, matrix + mm*k + (i/p)*nm, m, block);
					memcpy(matrix + mm*k + (i/p)*nm, block, size_m);
				}
				for(int k = 0; k < whole; k++)
				{
					matrixMulti(inv_block, m, m, id + mm*k + (i/p)*nm, m, block);
					memcpy(id + mm*k + (i/p)*nm, block, size_m);
				}
				memcpy(row_a, matrix + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
				memcpy(row_b, id + (i/p)*nm, __SIZEOF_DOUBLE__*nm);
			}
			else
			{
				for(int k = i + 1; k < whole; k++)
				{
					matrixMulti(inv_block, m, m, row_a + mm*k, m, block);
					memcpy(row_a + mm*k, block, size_m);
				}
				for(int k = 0; k < whole; k++)
				{
					matrixMulti(inv_block, m, m, row_b + mm*k, m, block);
					memcpy(row_b + mm*k, block, size_m);
				}
			}
			counter = 0;
			for(int l = my_rank; l < whole; l += p)
			{
				if(l != i)
				{
					for(int k = i + 1; k < whole; k++)
					{
						matrixMulti(matrix + (counter*nm + i*mm), m, m, row_a + k*mm, m, block);
						blockSubsruction(matrix + (counter*nm + k*mm), block, m, m);
					}
					for(int k = 0; k < whole; k++)
					{
						matrixMulti(matrix + (counter*nm + i*mm), m, m, row_b + k*mm, m, block);
						blockSubsruction(id + (counter*nm + k*mm), block, m, m);
					}
				}
				counter++;
			}
		}
	}
        for(int l = 0; l < whole; )
        {
                if(swap_in_matrix[l] != l)
                {
                        swapRowsInMatrixMPI(n, m, p, my_rank, swap_in_matrix[l], l, row_a, row_b, id);
                        swapRowsInArray(swap_in_matrix, swap_in_matrix[l], l);
                }
                else l++;
        }
	delete [] swap_in_block;
	delete [] swap_in_matrix;
	delete [] block;
	delete [] inv_block;
	return 0;
}
void makeRowsMatrix(double *matrix, double *e, int n, int m, int my_rank, int p)
{
    int whole = n/m;
    int rest = n % m;
    int current, l, count = 0, tag = 0;
    MPI_Status status;
    for(int i = 0; i < whole; i++)
	{
        current = my_rank;
        l = 0;
        while(current < whole)
        {
            if(i % p == my_rank)
            {
                memcpy(matrix + current*m*m + count*m*n, e + l*m*n + i*m*m,sizeof(double)*m*m);
                if(current + p - my_rank < whole)
                {
                    for(int j = 0; j < p; j++)
                    {
                        if(j != my_rank)
                        {
                            MPI_Recv(matrix + (j + p*l)*m*m + count*m*n, m*m, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                        }
                    }
                }
                else
                {
                    for(int j = 0; j < whole - current + my_rank; j++)
                    {

                        if(j != my_rank)
                        {
                            MPI_Recv(matrix + (j + p*l)*m*m + count*m*n, m*m, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                        }
                    }
                }
            }
            else
            {
                MPI_Send(e + l*m*n + i*m*m, m*m, MPI_DOUBLE, i % p, tag, MPI_COMM_WORLD);
            }
            l++;
            current += p;
        }
        if( i % p == my_rank)
        {
            for(int j = 0; j < whole - current + my_rank; j++)
            {
                if(j != my_rank)
                {
                    MPI_Recv(matrix + (j + p*l)*m*m + count*m*n, m*m, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                }
            }
        }
        if(rest != 0 && i % p == my_rank && whole % p == my_rank)
        {
            memcpy(matrix + whole*m*m + count*m*n, e + l*m*n + i*rest*m,sizeof(double)*rest*m);
        }
        else if(rest != 0 && i % p == my_rank)
        {
            MPI_Recv(matrix + whole*m*m + count*m*n, m*rest, MPI_DOUBLE, whole % p , tag, MPI_COMM_WORLD, &status);
        }
        else if(rest != 0 && whole % p == my_rank)
        {
            MPI_Send(e + l*m*n + i*rest*m, m*rest, MPI_DOUBLE, i % p, tag, MPI_COMM_WORLD);
        }
        if(i % p == my_rank)
        {
            count++;
        }
    }
    if(rest != 0)
    {
        current = my_rank;
        l = 0;
        while(current < whole)
        {
            if(whole % p == my_rank)
            {
                memcpy(matrix + current*m*rest + count*m*n, e + l*m*n + whole*m*m,sizeof(double)*m*rest);
                if(current + p - my_rank < whole)
                {
                    for(int j = 0; j < p; j++)
                    {
                        if(j != my_rank)
                        {
                            MPI_Recv(matrix + (j + p*l)*m*rest + count*m*n, m*rest, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                        }
                    }
                }
                else
                {
                    for(int j = 0; j < whole - current + my_rank; j++)
                    {

                        if(j != my_rank)
                        {
                            MPI_Recv(matrix + (j + p*l)*m*rest + count*m*n, m*rest, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                        }
                    }
                }
            }
            else
            {
                MPI_Send(e + l*m*n + whole*m*m, m*rest, MPI_DOUBLE, whole % p, tag, MPI_COMM_WORLD);

            }
            l++;
            current += p;
        }
        if( whole % p == my_rank)
        {
            for(int j = 0; j < whole - current + my_rank; j++)
            {
                if(j != my_rank)
                {
                    MPI_Recv(matrix + (j + p*l)*m*rest + count*m*n, m*rest, MPI_DOUBLE, j , tag, MPI_COMM_WORLD, &status);
                }
            }
        }
        if(rest != 0 && whole % p == my_rank && whole % p == my_rank)
        {
            memcpy(matrix + whole*m*rest + count*m*n, e + l*m*n + whole*rest*m,sizeof(double)*rest*rest);
        }
        else if(rest != 0 && whole % p == my_rank)
        {
            MPI_Recv(matrix + whole*m*rest + count*m*n, rest*rest, MPI_DOUBLE, whole % p , tag, MPI_COMM_WORLD, &status);
        }
        else if(rest != 0 && whole % p == my_rank)
        {
            MPI_Send(e + l*m*n + whole*rest*m, rest*rest, MPI_DOUBLE, whole % p, tag, MPI_COMM_WORLD);
        }
    }
}
void printMatrix(const double* matrix, int n, int m, int num_print)
{
        if(n < num_print)
        {
                cout<<"Печать матрицы невозможна, так как введены неверные параметры."<<endl;
                exit(1);
        }
        int remainder, old_row = 0, new_row = 0, pos = 0, amount = 1, counter = 0;
        remainder = n % m;
        while(counter < (num_print * num_print))
        {
                if(old_row % m == 0 && old_row > 0)
                {
                        new_row++;
                        pos = new_row * n * m;
                        amount = 1;
                        old_row=0;
                        cout<<setw(13)<<matrix[pos]<<" ";
                        counter++;
                        amount++;
                        pos++;
                        continue;
                }
                cout<<setw(13)<<matrix[pos]<<" ";
                counter++;
                if(amount == num_print)
                {
                        old_row++;
                        amount=1;
                        pos = n * m * new_row + m * old_row;
                        cout<<endl;
                        continue;
                }
                if(amount % m == 0)
                {
                        if((amount + remainder) != n || remainder == 0)
                        {
                                pos = n * m * new_row + amount * m + old_row * m;
                        }
                        else
                        {
                                pos = n * m * new_row + amount * m + remainder * old_row;
                        }
                        if(new_row * m + remainder == n)
                        {
                                pos = n * m * new_row + amount * remainder + old_row * m;
                        }
                        if(new_row * m + remainder == n && (amount + remainder) == n)
                        {
                                pos = n * m * new_row + amount * remainder + old_row * remainder;
                        }
                        amount++;
                        continue;
                }
                amount++;
                pos++;
        }
        cout<<endl;
        cout<<endl;
}
double residualMPI(double* matrix, double* id, double* row, int n, int m, int my_rank, int p, int max_rows, int file_or_formula, char* file_name)
{
    int whole = n/m;
    int remainder = n % m;
	int original_rank = my_rank;
	int variable_rank = my_rank;
	int mm = m*m;
	int nm = n*m;
	int rm = remainder*m;
	int k = 0;
	int source = (my_rank == 0) ? p - 1: my_rank - 1;
	int recipient = (my_rank == p - 1) ? 0: my_rank + 1;
    MPI_Status status;
    double *mult, *sum, *buf, *temp;
    mult = new double[m*m];
    sum = new double[m*m];
	buf = new double[n];
    int tag = 0;
    int counter_1 = 0;
    int counter_2 = 0;
    double max = 0.0;
	memset(row, 0, __SIZEOF_DOUBLE__*m*n);
	makeRowsMatrix(matrix, id, n, m, my_rank, p);
	createMatrix(id, row, n, m, p, my_rank, file_or_formula, file_name);
	memset(row, 0, __SIZEOF_DOUBLE__*m*n);
	temp = id;
	id = matrix;
	matrix = temp;
	for(int h = 0; h < p; h++)
	{
		counter_1 = 0;
		for(k = variable_rank; k < whole; k += p)
		{
			counter_2 = 0;
			for(int l = original_rank; l < whole; l += p)
			{
				memset(sum, 0, __SIZEOF_DOUBLE__*mm);
				for(int i = 0; i < whole; i++)
				{
					matrixMulti(matrix + i*mm + counter_1*nm, m, m, id + i*mm + counter_2*nm, m, mult);
					blockSum(sum, mult, m, m);
				}
				if(remainder != 0)
				{
					matrixMulti(matrix + whole*mm + counter_1*nm , m, remainder, id + whole*mm + counter_2*nm, m, mult);
					blockSum(sum, mult, m, m);
				}
				if(k == l) blockMinusId(sum, m);
				for(int j = 0; j < m; j++)
				{
					for(int i = 0; i < m; i++)
					{
						row[j + k*m] += fabs(sum[j*m + i]);
					}
				}
				counter_2++;
			}
			if(remainder != 0 && whole % p == original_rank)
			{
				memset(sum, 0, __SIZEOF_DOUBLE__*mm);
				for(int i = 0; i < whole; i++)
				{
					matrixMulti(matrix + i*mm + counter_1*nm, m, m, id + i*rm + counter_2*nm, remainder, mult);
					blockSum(sum, mult, m, remainder);
				}
				matrixMulti(matrix + whole*mm + counter_1*nm , m, remainder, id + whole*rm + counter_2*nm, remainder, mult);
			
				blockSum(sum, mult, m, remainder);
				for(int j = 0; j < m; j++)
				{
					for(int i = 0; i < remainder; i++)
					{
						row[j + k*m] += fabs(sum[j*remainder + i]);
					}
				}
			}
			counter_1++;
		}
		if(remainder != 0 && whole % p == variable_rank)
		{
			counter_2 = 0;
			for(int l = original_rank; l < whole; l += p)
			{
				memset(sum, 0, __SIZEOF_DOUBLE__*mm);
				for(int i = 0; i < whole; i++)
				{
					matrixMulti(matrix + i*rm + counter_1*nm, remainder, m, id + i*mm + counter_2*nm, m, mult);
					blockSum(sum, mult, remainder, m);
				}
				matrixMulti(matrix + whole*rm + counter_1*nm , remainder, remainder, id + whole*mm + counter_2*nm, m, mult);
				blockSum(sum, mult, remainder, m);
				for(int j = 0; j < remainder; j++)
				{
					for(int i = 0; i < m; i++)
					{
						row[j + whole*m] += fabs(sum[j*m + i]);
					}
				}
				counter_2++;
			}
			if(whole % p == original_rank)
			{
				memset(sum, 0, __SIZEOF_DOUBLE__*mm);
				for(int i = 0; i < whole; i++)
				{
					matrixMulti(matrix + i*rm + counter_1*nm, remainder, m, id + i*rm + counter_2*nm, remainder, mult);
					blockSum(sum, mult, remainder, remainder);
				}
				matrixMulti(matrix + whole*rm + counter_1*nm , remainder, remainder, id + whole*rm + counter_2*nm, remainder, mult);
				blockSum(sum, mult, remainder, remainder);
				blockMinusId(sum, remainder);
				for(int j = 0; j < remainder; j++)
				{
					for(int i = 0; i < remainder; i++)
					{
						row[j + whole*m] += fabs(sum[j*remainder + i]);
					}
				}
			}
		}
		MPI_Sendrecv_replace(matrix, nm*max_rows, MPI_DOUBLE, recipient, tag, source, tag, MPI_COMM_WORLD, &status);
		variable_rank = (variable_rank == 0) ? p - 1: variable_rank - 1;
	}
	if(original_rank == 0)
	{
		for(int i = 1; i < p; i++)
		{
			MPI_Recv(buf, n, MPI_DOUBLE, MPI_ANY_SOURCE , tag, MPI_COMM_WORLD, &status);
			for(int j = 0; j < n; j++)
			{
				row[j] += buf[j];
			}
		}
		max = row[0];
		for(int i = 1; i < n; i++) 
		{
			max = (max > row[i]) ? max: row[i];
		}
	}
	else MPI_Send(row, n, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	delete [] sum;
	delete [] mult;
	delete [] buf;
	return max;
}
