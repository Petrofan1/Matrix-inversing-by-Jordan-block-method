#include "Jordan_Header.hpp"
using namespace std;
void getBlock(const double* matrix, double* block, int g, int h, int n, int m)
{
        int k = (n/m) ;
        if(g == k)
        {
                if(h == k)
                {
                        memcpy(block, (matrix + ((g) * n * m + (h) * (n % m) * m)), sizeof(double) * (n % m) * (n % m));
                }
                else
                {
                        memcpy(block, (matrix + ((g) * n * m + (h) * (n % m) * m)), sizeof(double) * (n % m) * m);
                }
        }
        else
        {
                if( h == k)
                {
                        memcpy(block, (matrix + ((g) * n * m + (h) * m * m)), sizeof(double) * (n % m) * m);
                }
                else
                {
                        memcpy(block, (matrix + ((g) * n * m + (h) * m * m)), sizeof(double) * m * m);
                }
        }
}
void putBlock(double* matrix, const double* block, int g, int h, int n, int m)
{
        int k = (n/m);
        if(g == k)
        {
                if(h == k)
                {
                        memcpy((matrix + ((g) * n * m + (h) * (n % m) * m)), block,  sizeof(double) * (n % m) * (n % m));
                }
                else
                {
                        memcpy((matrix + ((g) * n * m + (h) * (n % m) * m)), block,  sizeof(double) * (n % m) * m);
                }
        }
        else
        {
                if(h == k)
                {
                        memcpy( (matrix + ((g) * n * m + (h) * m * m)), block, sizeof(double) * (n % m) * m);
                }
                else
                {
                        memcpy((matrix + ((g) * n * m + (h) * m * m)), block,  sizeof(double) * m * m);
                }
        }
}
double blockNorm(double* matrix, int rows, int columns)
{
        double max = 0;
        double current;
        for(int i = 0; i < rows; i++)
        {
                current = 0;
                for(int j = 0; j < columns; j++)
                {
                        current += fabs(matrix[i * columns + j]);
                }

                if(current > max)
                {
                        max = current;
                }
        }
        return max;
}
inline double matrixValue (int i, int j)
{
    return fabs(i-j);
    //return 1/double(i+j-1);
}
int readMatrixFromFile(double* matrix, int n, int m, char* file_name)
{
        FILE *file;
        file = fopen(file_name.c_str(), "r");       //Предполагается, что матрица в файле ровно того размера,
        if (file == NULL)                           //что был указан при запуске программы. Этот факт дополнительно не проверяется.
        {
                cout<<"Файл не найден."<<endl;
                exit(1);
        }
        double current;
        int counter = 0;
        int remainder, old_row = 0, new_row = 0, pos = 0, amount = 1;
        remainder = n % m;
        while(fscanf(file, "%lf", &current) == 1)
        {
                if(old_row % m == 0 && old_row > 0)
                {
                        new_row++;
                        pos = new_row * n * m;
                        amount = 1;
                        old_row=0;
                        matrix[pos] = current;
            counter++;
                        amount++;
                        pos++;
                        continue;
                }
                matrix[pos] = current;
        counter ++;
                if(amount == n)
                {
                        old_row++;
                        amount=1;
                        pos = n * m * new_row + m * old_row;
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
        fclose(file);
        if(counter != n*n)
        {
            cout<<"\nНе хватает данных"<<endl;
            return -1;
        }
        return 0;
}
void createIdMatrix(double* matrix, int n, int m)
{
        int remainder, old_row = 0, new_row = 0, pos = 0, amount = 1, counter = 0;
        remainder = n % m;
        while(counter < (n * n))
        {
                if(old_row % m == 0 && old_row > 0)
                {
                        new_row++;
                        pos = new_row * n * m;
                        amount = 1;
                        old_row = 0;
                        matrix[pos] = (new_row*m + 1 == amount) ? 1 : 0;
                        counter++;
                        amount++;
                        pos++;
                        continue;
                }
                matrix[pos] = ((new_row*m + old_row) + 1 == amount) ? 1 : 0;
                counter++;
                if(amount == n)
                {
                        old_row++;
                        amount=1;
                        pos = n * m * new_row + m * old_row;
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
}
void createMatrix(double* matrix, int n, int m)
{
        int remainder, old_row = 0, new_row = 0, pos = 0, amount = 1, counter = 0;
        remainder = n % m;
        while(counter < (n * n))
        {
                if(old_row % m == 0 && old_row > 0)
                {
                        new_row++;
                        pos = new_row * n * m;
                        amount = 1;
                        old_row = 0;
                        matrix[pos] = matrixValue(old_row + new_row * m + 1, amount);
                        counter++;
                        amount++;
                        pos++;
                        continue;
                }
                matrix[pos] = matrixValue(old_row + new_row * m + 1, amount);
                counter++;
                if(amount == n)
                {
                        old_row++;
                        amount=1;
                        pos = n * m * new_row + m * old_row;
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
void matrixMulti(double* A, int a1, int a2, double* B, int b2, double* Res)
{
        double *pB=B, *pA=A, *pRes=Res, sum[9];
        int remainder_i = a1 % 3;
        int remainder_j = b2 % 3;
        double a0, a9, a3, b0, b1, b3;
        for(int i = 0; i < a1 - remainder_i; i+=3)
        {
                pRes=Res + i*b2;
                for(int j = 0; j < b2 - remainder_j; j+=3)
                {
                    pA=A + i*a2;
                    pB=B  + j;
                        memset(sum, 0, sizeof(double)*9);
                        for(int k = 0; k < a2; k++)
                        {
                                a0=pA[0];
                                a9=pA[a2];
                                a3=pA[2*a2];
                                b0=pB[0];
                                b1=pB[1];
                                b3=pB[2];
                                sum[0]+=a0*b0;
                                sum[1]+=a0*b1;
                                sum[2]+=a0*b3;
                                sum[3]+=a9*b0;
                                sum[4]+=a9*b1;
                                sum[5]+=a9*b3;
                                sum[6]+=a3*b0;
                                sum[7]+=a3*b1;
                                sum[8]+=a3*b3;
                                pA++;
                                pB+=b2;
                        }
                        pRes[0]=sum[0];
                        pRes[1]=sum[1];
                        pRes[2]=sum[2];
                        pRes[b2]=sum[3];
                        pRes[b2+1]=sum[4];
                        pRes[b2+2]=sum[5];
                        pRes[2*b2]=sum[6];
                        pRes[2*b2+1]=sum[7];
                        pRes[2*b2+2]=sum[8];
                        pRes+=3;
                }
        }
        for(int i = a1 - remainder_i; i < a1; i++)
        {
                for(int j = b2 - remainder_j; j < b2; j++)
                {
                        sum[0]=0.;
                        pRes=Res + i*b2 + j;
                        for (int k = 0; k < a2; k++)
                        {
                                pA=A + i*a2 + k;
                                pB=B + k*b2 + j;
                                sum[0]+=pA[0]*pB[0];
                        }
                        pRes[0]=sum[0];
                }
        }
        for(int i = 0; i < a1 - remainder_i; i++)
        {
                for(int j = b2 - remainder_j; j < b2; j++)
                {
                        sum[0]=0.;
                        pRes=Res + i*b2 + j;
                        for (int k = 0; k < a2; k++)
                        {
                                pA=A + i*a2 + k;
                                pB=B + k*b2 + j;
                                sum[0]+=pA[0]*pB[0];
                        }
                        pRes[0]=sum[0];
                }
        }
        for(int i = a1 - remainder_i; i < a1; i++)
        {
                for(int j = 0; j < b2 - remainder_j; j++)
                {
                        sum[0]=0.;
                        pRes=Res + i*b2 + j;
                        for(int k = 0; k < a2; k++)
                        {
                                pA=A + i*a2 + k;
                                pB=B + k*b2 + j;
                                sum[0]+=pA[0]*pB[0];
                        }
                        pRes[0]=sum[0];
                }
        }
}
double getIJ(double* matrix, int i, int j, int n, int m)
{
        int k = n/m;
        int I = i/m;
        int J = j/m;
        if(J != k && I!=k)
        {
                return matrix[n*m*I+m*m*J+m*(i%m)+(j%m)];
        }
        if(J == k && I != k)
        {
                return matrix[n*m*I+m*m*J+(n%m)*(i%m)+(j%m)];
        }
        if(J != k && I == k)
        {
                return matrix[n*m*I+m*(n%m)*J+(n%m)*(i%m)+(j%m)];
        }
        return matrix[n*m*k+m*(n%m)*k+(n%m)*(i%m)+(j%m)];
}
void printBlock(double *block, int n, int m)
{
        for(int i = 0; i < n; i++)
        {
                for(int j = 0; j < m; j++)
                {
                        cout<<setw(12)<<block[j + i*m]<<" ";
                }
                cout<<endl;
        }
        cout<<endl;
}
void blockMulti(double* A, double* B, double* Res, int n, int m) //Умножение предусмотрено только для квадратных матриц одного размера
{
        int whole = n/m;
        int remainder = n%m;
        double *a, *b, *sum, *temp;
        a = new double[m*m];
        temp = new double[m*m];
        b = new double[m*m];
        sum = new double[m*m];
        if(remainder != 0)
        {
                for(int i = 0; i < whole; i++)
                {
                        for(int j = 0; j < whole; j++)
                        {
                                memset(sum, 0, m*m*sizeof(double));
                                for(int s = 0; s < whole; s++)
                                {
                                        memset(temp, 0, m*m*sizeof(double));
                                        getBlock(A, a, i, s, n, m);
                                        getBlock(B, b, s, j, n, m);
                                        matrixMulti(a, m, m, b, m, temp);
                                        for(int l = 0; l < m*m; l++)
                                        {
                                                sum[l]+=temp[l];
                                        }
                                }
                                memset(temp, 0, m*m*sizeof(double));
                                getBlock(A, a, i, whole, n, m);
                                getBlock(B, b, whole, j, n, m);
                                matrixMulti(a, m, remainder, b,  m, temp);
                                for(int l = 0; l < m*m; l++)
                                {
                                        sum[l]+=temp[l];
                                }
                                putBlock(Res, sum, i, j, n, m);
                        }
                }
                for(int j = 0; j < whole; j++)
                {
                        memset(sum, 0, m*m*sizeof(double));
                        for(int s = 0; s < whole; s++)
                        {
                                memset(temp, 0, m*m*sizeof(double));
                                getBlock(A, a, whole, s, n, m);
                                getBlock(B, b, s, j, n, m);
                                matrixMulti(a, remainder, m, b, m, temp);
                                for(int l = 0; l < (remainder)*m; l++)
                                {
                                        sum[l]+=temp[l];
                                }
                        }
                        memset(temp, 0, m*m*sizeof(double));
                        getBlock(A, a, whole, whole, n, m);
                        getBlock(B, b, whole, j, n, m);
                        matrixMulti(a, remainder, remainder, b, m, temp);
                        for(int l = 0; l < (remainder)*m; l++)
                        {
                                sum[l]+=temp[l];
                        }
                        putBlock(Res, sum, whole, j, n, m);
                }
                for(int i = 0; i < whole; i++)
                {
                        memset(sum, 0, m*m*sizeof(double));
                        for(int s = 0; s < whole; s++)
                        {
                                memset(temp, 0, m*m*sizeof(double));
                                getBlock(A, a, i, s, n, m);
                                getBlock(B, b, s, whole, n, m);
                                matrixMulti(a, m, m, b, remainder, temp);
                                for(int l = 0; l < (remainder)*m; l++)
                                {
                                        sum[l]+=temp[l];
                                }
                        }
                        memset(temp, 0, m*m*sizeof(double));
                        getBlock(A, a, i, whole, n, m);
                        getBlock(B, b, whole, whole, n, m);
                        matrixMulti(a, m, remainder, b,  remainder, temp);
                        for(int l = 0; l < (remainder)*m; l++)
                        {
                                sum[l]+=temp[l];
                        }
                        putBlock(Res, sum, i, whole, n, m);
                }
                memset(sum, 0, m*m*sizeof(double));
                for(int s = 0; s < whole; s++)
                {
                        memset(temp, 0, m*m*sizeof(double));
                        getBlock(A, a, whole, s, n, m);
                        getBlock(B, b, s, whole, n, m);
                        matrixMulti(a, remainder, m, b, remainder, temp);
                        for(int l = 0; l < (remainder)*(remainder); l++)
                        {
                                sum[l]+=temp[l];
                        }
                }
                memset(temp, 0, m*m*sizeof(double));
                getBlock(A, a, whole, whole, n, m);
                getBlock(B, b, whole, whole, n, m);
                matrixMulti(a, remainder, remainder, b, remainder, temp);
                for(int l = 0; l < (remainder)*(remainder); l++)
                {
                        sum[l]+=temp[l];
                }
                putBlock(Res, sum, whole, whole, n, m);
        }
        else
        {
                for(int i = 0; i < whole; i++)
                {
                        for(int j = 0; j < whole; j++)
                        {
                                memset(sum, 0, m*m*sizeof(double));
                                for(int s = 0; s < whole; s++)
                                {
                                        memset(temp, 0, m*m*sizeof(double));
                                        getBlock(A, a, i, s, n, m);
                                        getBlock(B, b, s, j, n, m);
                                        matrixMulti(a, m, m, b, m, temp);
                                        for(int l = 0; l < m*m; l++)
                                        {
                                                sum[l]+=temp[l];
                                        }
                                }
                                putBlock(Res, sum, i, j, n, m);
                        }
                }
        }
        delete [] a;
        delete [] b;
        delete [] sum;
        delete [] temp;
}
void createBlock(double *block, int n, int m)
{

        for(int i = 0; i < n; i++)
        {
                for(int j = 0; j < m; j++)
                {
                        block[j+i*m]=matrixValue(i,j);
                }
        }
}
void blockSubsruction(double* A, double* B, int n, int m)
{
        for(int i = 0; i < n; i++)
        {
                for(int j = 0; j < m; j++)
                {
                        A[i*m+j] -= B[i*m+j];
                }
        }
}
void swapColumn(double* block, int i, int j, int n)
{
        double temp;
        for(int k = 0; k < n; k++)
        {
                temp=block[k*n + j];
                block[k*n + j]=block[k*n + i];
                block[k*n + i]=temp;
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
void multiRow(double* block, double c, int i, int j, int n, int m)
{
        for(int k = j; k < n; k++)
        {
                block[i*m+k] *= c;
        }
}
void matrixId(double* block, int n)
{
        for(int s = 0; s < n; s++)
        {
                for(int p = 0; p < n; p++)
                {
                        block[s*n+p] = (s==p) ? 1 : 0;
                }
        }
}
void swapMatrixColumn(double* matrix, int i, int j, int n, int m, int whole, int remainder)
{
        double* temp;
        temp = new double[m*m];
        size_t size_m = m*m*(__SIZEOF_DOUBLE__);
        size_t size_r = m*remainder*(__SIZEOF_DOUBLE__);
        for (int k = 0;  k < whole; k++)
        {
                memcpy(temp, matrix + (k*n*m + j*m*m), size_m);
                memcpy(matrix + (k*n*m + j*m*m), matrix + (k*n*m + i*m*m), size_m);
                memcpy(matrix + (k*n*m + i*m*m), temp,  size_m);
        }
        if(remainder != 0)
        {
                memcpy(temp, matrix + (whole*n*m + j*m*remainder), size_r);
                memcpy(matrix + (whole*n*m + j*remainder*m), matrix + (whole*n*m + i*remainder*m), size_r);
                memcpy(matrix + (whole*n*m + remainder*i*m), temp,  size_r);
        }
        delete [] temp;
}
void swapRowsInArray(int* temp2, int i, int j)
{
        int s;
        s = temp2[i];
        temp2[i]=temp2[j];
        temp2[j] = s;
}
void swapMatrixRows(double* matrix, int i, int j, int n, int m, int whole, int remainder)
{
        double *temp;
        temp = new double[m*m];
        size_t size_m = m*m*(__SIZEOF_DOUBLE__);
        size_t size_r = m*remainder*(__SIZEOF_DOUBLE__);
        for (int k = 0;  k < whole; k++)
        {
                memcpy(temp, matrix + (i*n*m + k*m*m), size_m);
                memcpy(matrix + (i*n*m + k*m*m), matrix + (j*n*m + k*m*m), size_m);
                memcpy(matrix + (j*n*m + k*m*m), temp,  size_m);
        }
        if(remainder != 0)
        {
                memcpy(temp, matrix + (i*n*m + whole*m*m), size_r);
                memcpy(matrix + (i*n*m + whole*m*m), matrix + (j*n*m + whole*m*m), size_r);
                memcpy(matrix + (j*n*m + whole*m*m), temp,  size_r);
        }
        delete [] temp;
}
int blockInverse(double* block, int* temp2,  double* id, int n, double norm)
{

        int max;
        int s;
        double MAX;
        for(int i = 0; i < n; i++)
        {
                temp2[i]=i;
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
                                //cout<<MAX<<endl;
                        }
                }
                if(MAX < EPS*norm)
                {
                        return -1;
                }
                swapColumn(block, j, max, n);
                s = temp2[j];
                temp2[j]=temp2[max];
                temp2[max] = s;
                for(int l = 0; l < j; l++)
                {
                        for(int k = j+1; k < n; k++)
                        {
                                block[l*n+k] -= block[j*n+k]*(1/(block[j*n+j]))*block[l*n+j];
                        }
                        for(int k = 0; k < n; k++)
                        {
                                id[l*n+k] -= id[j*n+k]*(1/(block[j*n+j]))*block[l*n+j];
                        }
                }
                for(int l = j + 1; l < n; l++)
                {
                        for(int k = j+1; k < n; k++)
                        {
                                block[l*n+k] -= block[j*n+k]*(1/(block[j*n+j]))*block[l*n+j];
                        }
                        for(int k = 0; k < n; k++)
                        {
                                id[l*n+k] -= id[j*n+k]*(1/(block[j*n+j]))*block[l*n+j];
                        }
                }
                for(int k = j+1; k < n; k++)
                {
                        block[j*n+k] *= (1/(block[j*n+j]));
                }
                for(int k = 0; k < n; k++)
                {
                        id[j*n+k] *=(1/(block[j*n+j]));
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
double matrixNorm(double* matrix, int n, int m)
{
        double max = 0;
        double current;
        for(int i = 0; i < n; i++)
        {
                current = 0;
                for(int j = 0; j < n; j++)
                {
                        current += fabs(getIJ(matrix, i, j, n, m));
                }
                if(current > max)
                {
                        max = current;
                }
        }
        return max;
}
int matrixInverse(double* matrix, double* res, int n, int m, int* temp, int* temp2, double* block, double*  help1, double* inv_block)
{
        double matrix_norm = matrixNorm(matrix, n, m);
        int s;
        int whole = n/m;
        int nm = n*m;
        int mm = m*m;
        int remainder = n%m;
        int rm = remainder*m;
        int min;
        double MIN = __DBL_MAX__;
        double block_norm;
        size_t size_m = m*m*(__SIZEOF_DOUBLE__);
        if(remainder != 0)
        {
                for(int i = 0; i < whole; i++)
                {
                        min = i;
                        MIN = __DBL_MAX__;
                        for(int j = i; j < whole; j++)
                        {
                                getBlock(matrix, block, i, j, n, m);
                                if(blockInverse(block, temp, inv_block, m, matrix_norm) == -1)
                                {
                                        continue;
                                }
                                block_norm = blockNorm(inv_block, m, m);
                                if(block_norm < MIN)
                                {
                                        min = j;
                                        MIN = block_norm;
                                }
                        }
                        if((fabs(MIN - __DBL_MAX__))<EPS*matrix_norm)
                        {
                                cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
                                return -1;
                        }
                                swapMatrixColumn(matrix, i, min, n, m, whole, remainder);
                                s = temp2[i];
                                temp2[i]=temp2[min];
                                temp2[min] = s;
                        getBlock(matrix, block, i, i, n, m);
                        blockInverse(block, temp, inv_block, m, matrix_norm);
                        for(int k = i+1; k < whole; k++)
                        {
                                //memset(help1, 0, size_m);
                                matrixMulti(inv_block, m, m, matrix + (i*nm + mm*k),m, help1);
                                memcpy((matrix + ((i) * n * m + (k) * m * m)), help1,  size_m);
                        }
                        //memset(help1, 0, size_m);
                        matrixMulti(inv_block, m, m, matrix + (i*nm + whole*mm), remainder, help1);
                        putBlock(matrix, help1, i, whole, n, m);
                        for(int k = 0; k < whole; k++)
                        {
                                //memset(help1, 0, size_m);
                                matrixMulti(inv_block, m, m, res + (i*nm + mm*k), m, help1);
                                memcpy((res + ((i) * n * m + (k) * m * m)), help1,  size_m);
                        }
                        //memset(help1, 0, size_m);
                        matrixMulti(inv_block, m, m, res + (i*nm + whole*mm), remainder, help1);
                        putBlock(res, help1, i, whole, n, m);
                        for(int l = 0; l < i; l++)
                        {
                                for(int k = i+1; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + k*mm), m, block);
                                        blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
                                }
                               // memset(block, 0, size_m);
                                matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + whole*mm), remainder, block);
                                blockSubsruction(matrix + (l*nm + whole*mm), block, m, remainder);
                                for(int k = 0; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + k*mm), m, block);
                                        blockSubsruction(res + (l*nm + k*mm), block, m, m);
                                }
                               // memset(block, 0, size_m);
                                matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + whole*mm), remainder, block);
                                blockSubsruction(res + (l*nm + whole*mm), block, m, remainder);
                        }
                        for(int l = i + 1; l < whole; l++)
                        {
                                for(int k = i+1; k < whole; k++)
                                {
                                       // memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + k*mm), m, block);
                                        blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
                                }
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + whole*mm), remainder, block);
                                blockSubsruction(matrix + (l*nm + whole*mm), block, m, remainder);
                                for(int k = 0; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + k*mm), m, block);
                                        blockSubsruction(res + (l*nm + k*mm), block, m, m);
                                }
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + whole*mm), remainder, block);
                                blockSubsruction(res + (l*nm + whole*mm), block, m, remainder);
                        }
                        for(int k = i+1; k < whole; k++)
                        {
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (whole*nm + i*rm), remainder, m, matrix + (i*nm + k*mm),  m, block);
                                blockSubsruction(matrix + (whole*nm + k*rm), block, remainder, m);
                        }
                        //memset(block, 0, size_m);
                        matrixMulti(matrix + (whole*nm + i*rm), remainder, m, matrix + (i*nm + whole*mm), remainder, block);
                        blockSubsruction(matrix + (whole*nm + whole*rm), block, remainder, remainder);
                        for(int k = 0; k < whole; k++)
                        {
                               // memset(block, 0, size_m);
                                matrixMulti(matrix + (whole*nm + i*rm), remainder, m, res + (i*nm + k*mm), m, block);
                                blockSubsruction(res + (whole*nm + k*rm), block, remainder, m);
                        }
                        //memset(block, 0, size_m);
                        matrixMulti(matrix + (whole*nm + i*rm), remainder, m, res + (i*nm + whole*mm), remainder, block);
                        blockSubsruction(res + (whole*nm + whole*rm), block, remainder, remainder);
                }
                memcpy(block, (matrix + ((whole) * n * m + (whole) * (n % m) * m)), sizeof(double) * (n % m) * (n % m));
                if(blockInverse(block, temp, inv_block, remainder, matrix_norm) == -1)
                {
                        cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
                        return -1;
                }
                for(int k = 0; k < whole; k++)
                {
                        //memset(help1, 0, size_m);
                        matrixMulti(inv_block, remainder, remainder, res + (whole*nm + rm*k), m, help1);
                        putBlock(res, help1, whole, k, n, m);
                }
                //memset(help1, 0, size_m);
                matrixMulti(inv_block, remainder, remainder, res + (whole*nm + whole*rm),  remainder, help1);
                putBlock(res, help1, whole, whole, n, m);
                for(int l = 0; l < whole; l++)
                {
                        for(int k = 0; k < whole; k++)
                        {
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (l*nm + whole*mm), m, remainder, res + (whole*nm + k*rm),  m, block);
                                blockSubsruction(res + (l*nm + k*mm), block, m, m);
                        }
                        //memset(block, 0, size_m);
                        matrixMulti(matrix + (l*nm + whole*mm), m, remainder, res + (whole*nm + whole*rm), remainder, block);
                        blockSubsruction(res + (l*nm + whole*mm), block, m, remainder);
                }
        }
        else
        {
                for(int i = 0; i < whole; i++)
                {
                        min = i;
                        MIN = __DBL_MAX__;
                        for(int j = i; j < whole; j++)
                        {
                                getBlock(matrix, block, i, j, n, m);
                                if(blockInverse(block, temp, inv_block, m, matrix_norm) == -1)
                                {
                                        continue;
                                }
                                block_norm = blockNorm(inv_block, m, m);
                                if(block_norm < MIN)
                                {
                                        min = j;
                                        MIN = block_norm;
                                }
                        }
                        if((fabs(MIN - __DBL_MAX__))<EPS*matrix_norm)
                        {
                                cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
                                return -1;
                        }
                                swapMatrixColumn(matrix, i, min, n, m, whole, remainder);
                                s = temp2[i];
                                temp2[i]=temp2[min];
                                temp2[min] = s;
                        getBlock(matrix, block, i, i, n, m);
                        blockInverse(block, temp, inv_block, m, matrix_norm);
                        for(int k = i+1; k < whole; k++)
                        {
                                //memset(help1, 0, size_m);
                                matrixMulti(inv_block, m, m, matrix + (i*nm + mm*k), m, help1);
                                memcpy((matrix + ((i) * n * m + (k) * m * m)), help1,  size_m);
                        }
                        for(int k = 0; k < whole; k++)
                        {
                                //memset(help1, 0, size_m);
                                matrixMulti(inv_block, m, m, res + (i*nm + mm*k), m, help1);
                                memcpy((res + ((i) * n * m + (k) * m * m)), help1,  size_m);
                        }
                        for(int l = 0; l < i; l++)
                        {
                                for(int k = i+1; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + k*mm), m, block);
                                        blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
                                }
                                for(int k = 0; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + k*mm), m, block);
                                        blockSubsruction(res + (l*nm + k*mm), block, m, m);
                                }
                        }
                        for(int l = i + 1; l < whole; l++)
                        {
                                for(int k = i+1; k < whole; k++)
                                {
                                       // memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + (i*nm + k*mm), m, block);
                                        blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
                                }
                                for(int k = 0; k < whole; k++)
                                {
                                        //memset(block, 0, size_m);
                                        matrixMulti(matrix + (l*nm + i*mm), m, m, res + (i*nm + k*mm), m, block);
                                        blockSubsruction(res + (l*nm + k*mm), block, m, m);
                                }
                        }
                        for(int k = i+1; k < whole; k++)
                        {
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (whole*nm + i*rm), remainder, m, matrix + (i*nm + k*mm), m, block);
                                blockSubsruction(matrix + (whole*nm + k*rm), block, remainder, m);
                        }
                        for(int k = 0; k < whole; k++)
                        {
                                //memset(block, 0, size_m);
                                matrixMulti(matrix + (whole*nm + i*mm), remainder, m, res + (i*nm + k*mm), m, block);
                                blockSubsruction(res + (whole*nm + k*mm), block, remainder, m);
                        }
                }
        }
        for(int l = 0; l < whole; )
        {
                if(temp2[l] != l)
                {
                        swapMatrixRows(res, temp2[l], l, n, m, whole, remainder);
                        swapRowsInArray(temp2, temp2[l], l);
                }
                else l++;
        }
        return 0;
}
double discrepansy2(double* matrix, double* inv_matrix,  int n, int m)
{

        int whole = n/m;
        int remainder = n%m;
        int nm = n*m;
        int mm = m*m;
        int rm = remainder*m;
        double *A=matrix, *B=inv_matrix;
        double *r;
        double *p;
        double *id;
        double *id2;
        double* m3;
        p = new double[whole + 1];
        r = new double[m];
        id = new double[m*m];
        id2 = new double[remainder*remainder];
        m3 = new double[m*m];
        matrixId(id, m);
        matrixId(id2, remainder);
        memset(r, 0, sizeof(double)*m);
        memset(p, 0, sizeof(double)*(whole+1));
        int i,j;
        for(i = 0; i < whole; i++)
        {
                for(j = 0; j < whole; j++)
                {
                        memset(m3, 0, sizeof(double)*m*m);
                        for(int s = 0; s < whole; s++)
                        {
                                matrixMulti(A + (i*nm + s*mm), m, m, B + (s*nm + j*mm), m, m3);
                        }
                        matrixMulti(A + (i*nm + whole*mm), m, remainder, B + (whole*nm + j*rm), m, m3);
                        if(i==j) blockSubsruction(m3, id, m, m);
                        for(int k = 0; k < m; k++)
                        {
                                for(int l = 0; l < m; l++)
                                {
                                        r[k] += fabs(m3[k*m + l]);
                                }
                        }

                }
                memset(m3, 0, sizeof(double)*m*m);
                for(int s = 0; s < whole; s++)
                {
                        matrixMulti(A + (i*nm + s*mm), m, m, B + (s*nm + whole*mm), remainder, m3);
                }
                matrixMulti(A + (i*nm + whole*mm), m, remainder, B + (whole*nm + j*rm), remainder, m3);
                for(int k = 0; k < m; k++)
                {
                        for(int l = 0; l < remainder; l++)
                        {
                                r[k] += fabs(m3[k*remainder + l]);
                        }
                }
                for(int v = 0; v < m; v++)
                {
                        p[i] = (p[i] > r[v]) ? p[i]: r[v];

                }
                memset(r, 0, sizeof(double)*m);
        }
        for(int q = 0; q < whole; q++)
        {
                memset(m3, 0, sizeof(double)*m*m);
                for(int s = 0; s < whole; s++)
                {
                        matrixMulti(A + (whole*nm + s*rm), remainder, m, B + (s*nm + q*mm), m, m3);
                }
                matrixMulti(A + (whole*nm + whole*rm), remainder, remainder, B + (whole*nm + q*rm), m, m3);
                for(int k = 0; k < remainder; k++)
                {
                        for(int l = 0; l < m; l++)
                        {
                                r[k] += fabs(m3[k*m + l]);
                        }
                }
        }
        memset(m3, 0, sizeof(double)*m*m);
        for(int s = 0; s < whole; s++)
        {
                matrixMulti(A + (whole*nm + s*rm), remainder, m, B + (s*nm + whole*mm), remainder, m3);
        }
        matrixMulti(A + (whole*nm + whole*rm), remainder, remainder, B + (whole*nm + whole*rm), remainder, m3);
        blockSubsruction(m3, id2, remainder, remainder);
        for(int k = 0; k < remainder; k++)
        {
                for(int l = 0; l < remainder; l++)
                {
                        r[k] += fabs(m3[k*remainder + l]);
                }
        }
        for(int v = 0; v < remainder; v++)
        {
                p[whole] = (p[whole] > r[v]) ? p[whole]: r[v];
        }
        double max = 0.;
        for(int f = 0; f < whole + 1; f++)
        {
                max = (p[f] > max) ? p[f]: max;
        }
        delete [] r;
        delete [] p;
        delete [] id;
        delete [] id2;
        delete [] m3;
        return max;
}
double discrepansy(double *oldMatrix,int n,int m, double *matrix,double*row){
    int whole = int(n/m);
    int rest = n%m;
    double max = -1;
    int mm = m*m, nm = n*m, restm = rest*m,start1,start;
    double *c,*sum;
    c = new double[mm];
    sum = new double[mm];
    //double *row = new double[m];

    for(int k = 0; k < whole; k++){
        memset(row,0,sizeof (double) *m);
        for(int j = 0; j < whole;j++){

            memset(sum,0,sizeof(double)*mm);
            start = k*nm;
            start1 = j*mm;
            for(int i = 0; i< whole; i++){
                //getBlock(oldMatrix,start,m,m,a);
                //getBlock(matrix,start1,m,m,b);
                matrixMulti(oldMatrix + start,m,m,matrix + start1,m,c);
               // printM(c,m,m);
                for(int g = 0; g < mm; g++){
                    sum[g] += c[g];
                }
                start += mm;
                start1 += nm;
            }
            if(rest != 0){
                //getBlock(oldMatrix,start,m,rest,a);
                //getBlock(matrix,j*restm + whole*nm,rest,m,b);
                matrixMulti(oldMatrix + start,m,rest,matrix + j*restm + whole*nm,m,c);
                for(int l = 0; l < mm; l++){
                    sum[l] += c[l];
                }
            }
            if(k == j){
                for(int l = 0; l < m; l++){
                    sum[l*m + l] -=1;
                }
            }



          // putBlock(result,k*nm + j*mm,m,m,sum);
           for(int g =0; g < m; g++){
               for(int j =0; j < m; j++){
                   row[g] += fabs(sum[g*m + j]);
               }
           }
          if(rest != 0){
           memset(sum,0,sizeof(double)*mm);
           start = k*nm;
           start1 = whole*mm;
           for(int i = 0; i< whole; i++){
               //getBlock(oldMatrix,start,m,m,a);
               //getBlock(matrix,start1,m,rest,b);
               matrixMulti(oldMatrix + start,m,m,matrix + start1,rest,c);
               for(int l = 0; l < restm; l++){
                   sum[l] += c[l];
               }
               start += mm;
               start1 += nm;

           }

           //getBlock(oldMatrix,start,m,rest,a);
           //getBlock(matrix,whole*nm + whole*restm,rest,rest,b);
           matrixMulti(oldMatrix+start,m,rest,matrix + whole*nm + whole*restm,rest,c);
           for(int l = 0; l < restm; l++){
               sum[l] += c[l];
           }
          // putBlock(result,k*nm + whole*mm,m,m,sum);
           for(int g =0; g < m; g++){
               for(int j =0; j < rest; j++){
                   row[g] += fabs(sum[g*m + j]);
               }
           }

          }


        }


        for(int d =0; d < m; d++){
            if(max < row[d]){
                max = row[d];
            }
        }

    }

    if(rest != 0){
        memset(row,0,sizeof (double) *m);
        for(int k = 0; k < whole; k++){
            memset(sum,0,sizeof(double)*mm);
            start = whole*nm;
            start1 = k*mm;
            for(int i = 0; i< whole; i++){
                //getBlock(oldMatrix,start,rest,m,a);
                //getBlock(matrix,start1,m,m,b);
                matrixMulti(oldMatrix + start,rest,m,matrix + start1,m,c);
                for(int l = 0; l < restm; l++){
                    sum[l] += c[l];
                }
                start += restm;
                start1 += nm;
            }

            //getBlock(oldMatrix,start,rest,rest,a);
           // getBlock(matrix,k*restm + whole*nm,rest,m,b);
            matrixMulti(oldMatrix+start,rest,rest,matrix + k*restm + whole*nm,m,c);
            for(int l = 0; l < restm; l++){
                sum[l] += c[l];
            }
            //putBlock(result,whole*nm + k*restm,m,m,sum);
            for(int g =0; g < m; g++){
                for(int j =0; j < rest; j++){
                    row[g] += sum[g*m + j];
                }
            }

        }
        memset(sum,0,sizeof(double)*mm);

        for(int i = 0; i < whole; i++){
           // getBlock(oldMatrix,whole*nm + i*restm ,rest,m,a);
           // getBlock(matrix,whole*mm + i*nm,m,rest,b);
            matrixMulti(oldMatrix + whole*nm + i*restm,rest,m,matrix + whole*mm + i*nm,rest,c);
            for(int l = 0; l < rest*rest; l++){
                sum[l] += c[l];
            }
        }
        //getBlock(oldMatrix,whole*nm + whole*restm ,rest,rest,a);
        //getBlock(matrix,whole*nm + whole*restm,rest,rest,b);
       matrixMulti(oldMatrix + whole*nm + whole*restm,rest,rest,matrix + whole*nm + whole*restm,rest,c);
        for(int l = 0; l < rest*rest; l++){
            sum[l] += c[l];
        }
        for(int l = 0; l < rest; l++){
            sum[l*m + l] -=1;
        }
       // putBlock(result,whole*nm + whole*restm,rest,rest,sum);
        for(int g =0; g < rest; g++){
            for(int j =0; j < rest; j++){
                row[g] += sum[g*m + j];
            }
        }
        for(int d =0; d < rest; d++){
            if(max < row[d]){
                max = row[d];
            }
        }

    }
    //delete [] row;
    delete [] c;
    delete [] sum;
    return max;

}
