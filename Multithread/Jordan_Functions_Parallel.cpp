#include"Jordan_Header.hpp"
inline double matrixValue (int i, int j)
{
    return fabs(i-j);
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
double discrepansy(double *oldMatrix, double *matrix, int n,int m, double*row)
{
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
			memset(m3, 0, sizeof(double)*mm);
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
		memset(m3, 0, sizeof(double)*mm);
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
		memset(m3, 0, sizeof(double)*mm);
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
int readMatrixFromFile(double* matrix, int n, int m, string file_name)
{
	FILE *file;
	file = fopen(file_name.c_str(), "r");       //Предполагается, что матрица в файле ровно того размера,
	if (file == NULL)                           //что был указан при запуске программы. Этот факт дополнительно не проверяется.
	{
		cout<<"Cannot open file"<<endl;
		return -1;
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
		for(int k = j+1; k < n; k++)
		{
			block[j*n+k] /=(block[j*n+j]);
		}
		for(int k = 0; k < n; k++)
		{
			id[j*n+k] /=(block[j*n+j]);
		}
		for(int l = 0; l < j; l++)
		{
			for(int k = j+1; k < n; k++)
			{
				block[l*n+k] -= (block[j*n+k])*block[l*n+j];
			}
			for(int k = 0; k < n; k++)
			{
				id[l*n+k] -= (id[j*n+k])*block[l*n+j];
			}
		}
		for(int l = j + 1; l < n; l++)
		{
			for(int k = j+1; k < n; k++)
			{
				block[l*n+k] -= (block[j*n+k])*block[l*n+j];
			}
			for(int k = 0; k < n; k++)
			{
				id[l*n+k] -= (id[j*n+k])*block[l*n+j];
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
void blockSubsruction(double* A, double* B, int n, int m)
{
        for(int i = 0; i < n*m; i++)
        {
            A[i] -= B[i];
        }
}
void blockSum(double* A, double* B, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++)
		{
			A[i*m+j] += B[i*m+j];
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
		return matrix[n*m*I+m*(n%m)*J+(m)*(i%m)+(j%m)];
	}
	return matrix[n*m*k+m*(n%m)*k+(n%m)*(i%m)+(j%m)];
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
void matrixThreadsZeros(arg* a)
{
	int n = a->n;
	int m = a->m;
	int whole = n/m;
	int p = a->p;
	for(int i = a->num; i < whole; i+=p)
	{
		memset(a->a + n*m*i, 0, n*m*(__SIZEOF_DOUBLE__));
		memset(a->b + n*m*i, 0, n*m*(__SIZEOF_DOUBLE__));
	}
	if((whole - a->num) % p == 0 && n%m != 0)
	{
		memset(a->a + n*m*whole, 0, n*(n%m)*(__SIZEOF_DOUBLE__));
		memset(a->b + n*m*whole, 0, n*(n%m)*(__SIZEOF_DOUBLE__));
	}
	pthread_barrier_wait(a->barrier);
}
int matrixThreadsData(arg* a)
{
	static int error = 0;
    int m = a->m;
    int n = a->n;
    double *A = a->a;
    double *B = a->b;
    char* file_name = a->file_name;
	if(a->num == 0)
	{
		createIdMatrix(B, n, m);
		if(a->file_name == NULL)
		{
			createMatrix(A, n, m);
		}
		else
		{
			if(readMatrixFromFile(A, n, m, file_name) == -1) error = -2;
		}
	}
	pthread_barrier_wait(a->barrier);
	return error;
}
void matrixThreadsData2(arg* a)
{
    int m = a->m;
    int n = a->n;
    double *A = a->a;
    char* file_name = a->file_name;
	if(a->num == 0)
	{
		if(a->file_name == NULL)
		{
			createMatrix(A, n, m);
		}
		else
		{
			readMatrixFromFile(A, n, m, file_name);
		}
	}
	pthread_barrier_wait(a->barrier);
}
void matrixThreadsNorm(arg* a, double* norm)
{
    int m = a->m;
    int n = a->n;
    double *A = a->a;
    static double current_norm = 0.0;
    if(a->num == 0)
    {
        current_norm = matrixNorm(A, n, m);
    }
	pthread_barrier_wait(a->barrier);
    *norm = current_norm;
}
void findMainElement(arg* a, int i, double norm, int* index, double* min_norm, int* temp, double* inv_block, double* block)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    int p = a->p;
    int n = a->n;
    int m = a->m;
    int whole = n/m;
    double *A = a->a;
    static int t_in = 0;
    static int t_out = 0;
    static int current_index = 0;
    static double current_min_norm = __DBL_MAX__;
    int min = i;
    double MIN = __DBL_MAX__;
    double block_norm = 0;
    double start = a->num;
    while(start < i)
    {
            start += p;
    }
    pthread_mutex_lock(&mut);
    for(int j = start; j < whole; j+=p)
    {
            memcpy(block, A + i*n*m + j*m*m, m*m*__SIZEOF_DOUBLE__);
            if(blockInverse(block, temp, inv_block, m, norm) == -1)
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
    if(MIN < current_min_norm)
    {
            current_min_norm = MIN;
            current_index = min;
    }
    t_in++;
    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while(t_in < p)
        {
            pthread_cond_wait(&c_in, &mut);
        }
    }
        *min_norm = current_min_norm;
        *index = current_index;
        t_out++;
    if(t_out >= p)
    {
        current_min_norm = __DBL_MAX__;
        current_index = 0;
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while(t_out < p)
        {
            pthread_cond_wait(&c_out, &mut);
        }
    }
    pthread_mutex_unlock(&mut);
}
void swapMatrixColumn(arg* a, double* matrix, double* temp, int i, int j, int n, int m, int whole, int remainder, int p)
{
	size_t size_m = m*m*(__SIZEOF_DOUBLE__);
	size_t size_r = m*remainder*(__SIZEOF_DOUBLE__);
	for (int k = a->num;  k < whole; k+=p)
	{
		memcpy(temp, matrix + (k*n*m + j*m*m), size_m);
		memcpy(matrix + (k*n*m + j*m*m), matrix + (k*n*m + i*m*m), size_m);
		memcpy(matrix + (k*n*m + i*m*m), temp,  size_m);
	}
	if((whole - a->num) % p == 0 && remainder != 0)
	{
		memcpy(temp, matrix + (whole*n*m + j*m*remainder), size_r);
		memcpy(matrix + (whole*n*m + j*remainder*m), matrix + (whole*n*m + i*remainder*m), size_r);
		memcpy(matrix + (whole*n*m + remainder*i*m), temp,  size_r);
	}
}
void swapRowsInArray(int* temp2, int i, int j)
{
	int s;
	s = temp2[i];
	temp2[i]=temp2[j];
	temp2[j] = s;
}
void swapMatrixRows(double* temp, double* matrix, int i, int j, int n, int m, int whole, int remainder)
{
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
}
double getCPUTime()
{
    struct rusage t;
    getrusage(RUSAGE_THREAD, &t);
    return (double)t.ru_utime.tv_sec + (double)t.ru_utime.tv_usec/1000000.;
}
double getTime()
{
    struct timeval t;
    gettimeofday(&t, 0);
    return (double)t.tv_sec+(double)t.tv_usec/1000000.;
}
void matrixThreadsInverse(arg* a, double matrix_norm, double* block, double* help1, int* temp, int *temp2, double* inv_block)
{
        double MIN = __DBL_MAX__;
	int n = a->n;
	int m = a->m;
	int whole = n/m;
	int remainder = n%m;
	int nm = n*m;
	int mm = m*m;
	int rm = remainder*m;
	int p = a->p;
	size_t size_m = m*m*(__SIZEOF_DOUBLE__);
	int s;
        double* matrix = a->a;
        double* res = a->b;
	int min;
        double time = 0.;
	if(remainder != 0)
	{
		for(int i = 0; i < whole; i++)
		{
                        findMainElement(a, i, matrix_norm, &min, &MIN, temp, inv_block, block);
			if((fabs(MIN - __DBL_MAX__)) < EPS*matrix_norm)
			{
				cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
				a->error = -2;
			}
			if(min!=i)
			{
				swapMatrixColumn(a, a->a, block, i, min, n, m, whole, remainder, p);
				if(a->num == 0)
				{
					s = temp2[i];
					temp2[i] = temp2[min];
					temp2[min] = s;
				}
			}
                        pthread_barrier_wait(a->barrier);
                        getBlock(a->a, block, i, min, n, m);
                        blockInverse(block, temp, inv_block, m, matrix_norm);
                        int start = a->num;
                        while(start < i + 1)
                        {
                                start += p;
                        }
                        for(int k = start; k < whole; k+=p)
                        {
                                matrixMulti(inv_block, m, m, matrix + i*nm + mm*k,m, help1);
                                memcpy(matrix + i*nm + mm*k, help1,  size_m);
                        }
                        if((whole - a->num) % p == 0)
                        {
                            matrixMulti(inv_block, m, m, matrix + i*nm + whole*mm, remainder, help1);
                            putBlock(matrix + i*nm, help1, 0, whole, n, m);
                        }
                        for(int k = a->num; k < whole; k+=p)
                        {
                                matrixMulti(inv_block, m, m, res + i*nm + mm*k, m, help1);
                                memcpy(res + i*nm + (k) * m * m, help1,  size_m);
                        }
                        if((whole - a->num) % p == 0)
                        {
                                matrixMulti(inv_block, m, m, res + i*nm + whole*mm, remainder, help1);
                                putBlock(res + i*nm, help1, 0, whole, n, m);
                        }
                        pthread_barrier_wait(a->barrier);
                        for(int l = a->num; l < whole; l+=p)
			{
				if(l != i)
				{
					for(int k = i+1; k < whole; k++)
					{
						matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + i*nm + k*mm, m, block);
						blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
					}
					matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + i*nm + whole*mm, remainder, block);
					blockSubsruction(matrix + (l*nm + whole*mm), block, m, remainder);
					for(int k = 0; k < whole; k++)
					{
						matrixMulti(matrix + (l*nm + i*mm), m, m, res + i*nm + k*mm, m, block);
						blockSubsruction(res + (l*nm + k*mm), block, m, m);
					}
					matrixMulti(matrix + (l*nm + i*mm), m, m, res + i*nm + whole*mm, remainder, block);
					blockSubsruction(res + (l*nm + whole*mm), block, m, remainder);
				}
			}
			if((whole - a->num) % p == 0)
			{
				for(int k = i+1; k < whole; k++)
				{
					matrixMulti(matrix + (whole*nm + i*rm), remainder, m, matrix + i*nm + k*mm,  m, block);
					blockSubsruction(matrix + (whole*nm + k*rm), block, remainder, m);
				}
				matrixMulti(matrix + (whole*nm + i*rm), remainder, m, matrix + i*nm + whole*mm, remainder, block);
				blockSubsruction(matrix + (whole*nm + whole*rm), block, remainder, remainder);
				for(int k = 0; k < whole; k++)
				{
					matrixMulti(matrix + (whole*nm + i*rm), remainder, m, res + i*nm + k*mm, m, block);
					blockSubsruction(res + (whole*nm + k*rm), block, remainder, m);
				}
				matrixMulti(matrix + (whole*nm + i*rm), remainder, m, res + i*nm + whole*mm, remainder, block);
				blockSubsruction(res + (whole*nm + whole*rm), block, remainder, remainder);
			}
			pthread_barrier_wait(a->barrier);
		}
		if((whole - a->num) % p == 0)
		{
			memcpy(block, (matrix + ((whole) * n * m + (whole) * (n % m) * m)), sizeof(double) * (n % m) * (n % m));
			if(blockInverse(block, temp, inv_block, remainder, matrix_norm) == -1)
			{
				cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
				a-> error = -2;
			}
			for(int k = 0; k < whole; k++)
			{
				matrixMulti(inv_block, remainder, remainder, res + (whole*nm + rm*k), m, help1);
				putBlock(res, help1, whole, k, n, m);
			}
			matrixMulti(inv_block, remainder, remainder, res + (whole*nm + whole*rm),  remainder, help1);
			putBlock(res, help1, whole, whole, n, m);
		}
		pthread_barrier_wait(a->barrier);
		for(int l = a->num; l < whole; l+=p)
		{
			for(int k = 0; k < whole; k++)
			{
				matrixMulti(matrix + (l*nm + whole*mm), m, remainder, res + (whole*nm + k*rm),  m, block);
				blockSubsruction(res + (l*nm + k*mm), block, m, m);
			}
			matrixMulti(matrix + (l*nm + whole*mm), m, remainder, res + (whole*nm + whole*rm), remainder, block);
			blockSubsruction(res + (l*nm + whole*mm), block, m, remainder);
		}
		pthread_barrier_wait(a->barrier);
    }
	else
	{
		for(int i = 0; i < whole; i++)
                {
                    double time2 = getTime();
                        findMainElement(a, i, matrix_norm, &min, &MIN, temp, inv_block, block);
                        time += getTime() - time2;
			if((fabs(MIN - __DBL_MAX__)) < EPS*matrix_norm)
			{
				cout<<"\nМетод неприменим для заданного m или матрица необратима"<<endl;
				a->error = -2;
			}
			if(min!=i)
			{
				swapMatrixColumn(a, a->a, block, i, min, n, m, whole, remainder, p);
				if(a->num == 0)
				{	
					s = temp2[i];
					temp2[i] = temp2[min];
					temp2[min] = s;
				}
			}
                        pthread_barrier_wait(a->barrier);
                        getBlock(a->a, block, i, i, n, m);
			blockInverse(block, temp, inv_block, m, matrix_norm);
			int start = a->num;
			while(start < i + 1)
			{
				start += p;
			}
			for(int k = start; k < whole; k+=p)
			{
                                matrixMulti(inv_block, m, m, matrix + i*nm + mm*k, m, help1);
                                memcpy(matrix + i*nm + (k) * m * m, help1,  size_m);
			}
			for(int k = a->num; k < whole; k+=p)
			{
                                matrixMulti(inv_block, m, m, res + i*nm + mm*k, m, help1);
                                memcpy((res + i*nm + (k) * m * m), help1,  size_m);
			}
                        pthread_barrier_wait(a->barrier);
                        for(int l = a->num; l < whole; l+=p)
			{
				if(l != i)
				{
					for(int k = i+1; k < whole; k++)
					{
                                                matrixMulti(matrix + (l*nm + i*mm), m, m, matrix + i*nm  + k*mm, m, block);
						blockSubsruction(matrix + (l*nm + k*mm), block, m, m);
					}
					for(int k = 0; k < whole; k++)
					{
                                                matrixMulti(matrix + (l*nm + i*mm), m, m, res + i*nm + k*mm, m, block);
						blockSubsruction(res + (l*nm + k*mm), block, m, m);
					}
				}
			}
			pthread_barrier_wait(a->barrier);
		}
	}
	if(a->num == 0)
	{
		for(int l = 0; l < whole; )
		{
			if(temp2[l] != l)
			{
				swapMatrixRows(block, res, temp2[l], l, n, m, whole, remainder);
				swapRowsInArray(temp2, temp2[l], l);
			}
			else l++;
		}
	}
}
void max(arg* a, double max, double* MAXIMUM)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static double current_max = 0;
    int p = a->p;
    pthread_mutex_lock(&mut);
    {
        if(max > current_max)
        {
            current_max = max;
        }
    }
    t_in++;
        if(t_in >= p)
        {
            t_out = 0;
            pthread_cond_broadcast(&c_in);
        }
        else
        {
            while(t_in < p)
            {
                pthread_cond_wait(&c_in, &mut);
            }
        }
            *MAXIMUM = current_max;
            t_out++;
        if(t_out >= p)
        {
            t_in = 0;
            pthread_cond_broadcast(&c_out);
        }
        else
        {
            while(t_out < p)
            {
                pthread_cond_wait(&c_out, &mut);
            }
        }
        pthread_mutex_unlock(&mut);
}
void discrepansyThread(arg* a, double time)
{
	int n = a->n;
	int m = a->m;
	int whole = n/m;
	int p = a->p;
	int remainder = n%m;
	int nm = n*m;
	int mm = m*m;
	int rm = remainder*m;
        //double max = 0;
        double MAX = 0;
	double *A= a->a, *B=a->b;
	double *r;
	double *id;
	double *id2;
	double* m3;
	double* m4;
        double MAXIMUM;
	r = new double[m];
	id = new double[m*m];
	id2 = new double[remainder*remainder];
	m3 = new double[m*m];
	m4 = new double[m*m];
	matrixId(id, m);
	matrixId(id2, remainder);
	memset(r, 0, sizeof(double)*m);
	int i,j;
	for(i = a->num; i < whole; i+=p)
	{
		int counter = 1;
		j = i;
		while(counter != whole)
		{
			memset(m4, 0, sizeof(double)*mm);
			for(int s = 0; s < whole; s++)
			{
				matrixMulti(A + (i*nm + s*mm), m, m, B + (s*nm + j*mm), m, m3);
				blockSum(m4, m3, m, m);
			}
			if(remainder != 0)
			{
				matrixMulti(A + (i*nm + whole*mm), m, remainder, B + (whole*nm + j*rm), m, m3);
				blockSum(m4, m3, m, m); 
			}
			if(i==j) blockSubsruction(m4, id, m, m);
			for(int k = 0; k < m; k++)
			{
				for(int l = 0; l < m; l++)
				{
					r[k] += fabs(m4[k*m + l]);
				}
			}
			j++;
			if(j == whole) j = 0;
			counter++;
		}
		if(remainder != 0)
		{
			memset(m4, 0, sizeof(double)*mm);
			for(int s = 0; s < whole; s++)
			{
				matrixMulti(A + (i*nm + s*mm), m, m, B + (s*nm + whole*mm), remainder, m3);
				blockSum(m4, m3, m, remainder);
			}
			matrixMulti(A + (i*nm + whole*mm), m, remainder, B + (whole*nm + whole*rm), remainder, m3);
			blockSum(m4, m3, m, remainder);
			for(int k = 0; k < m; k++)
			{
				for(int l = 0; l < remainder; l++)
				{
					r[k] += fabs(m4[k*remainder + l]);
				}
			}
		}
		for(int v = 0; v < m; v++)
		{
			MAX = (MAX > r[v]) ? MAX: r[v];
		}
                memset(r, 0, sizeof(double)*m);
	}
	memset(m3, 0, m*m*__SIZEOF_DOUBLE__);
	if(remainder != 0 && (a->num - whole) % p == 0)
	{
		for(int q = 0; q < whole; q++)
		{
			memset(m4, 0, sizeof(double)*mm);
			for(int s = 0; s < whole; s++)
			{
				matrixMulti(A + (whole*nm + s*rm), remainder, m, B + (s*nm + q*mm), m, m3);
				blockSum(m4, m3, remainder, m);
			}
			matrixMulti(A + (whole*nm + whole*rm), remainder, remainder, B + (whole*nm + q*rm), m, m3);
			blockSum(m4, m3, remainder, m);
                        for(int k = 0; k < remainder; k++)
			{
				for(int l = 0; l < m; l++)
				{
					r[k] += fabs(m4[k*m + l]);
				}
			}
		}
		memset(m4, 0, sizeof(double)*m*m);
		for(int s = 0; s < whole; s++)
		{
			matrixMulti(A + (whole*nm + s*rm), remainder, m, B + (s*nm + whole*mm), remainder, m3);
			blockSum(m4, m3, remainder, remainder);
		}
		matrixMulti(A + (whole*nm + whole*rm), remainder, remainder, B + (whole*nm + whole*rm), remainder, m3);
		blockSum(m4, m3, remainder, remainder);
		blockSubsruction(m4, id2, remainder, remainder);
		for(int k = 0; k < remainder; k++)
		{
			for(int l = 0; l < remainder; l++)
			{
				r[k] += fabs(m4[k*remainder + l]);
			}
		}
		for(int v = 0; v < remainder; v++)
		{
			MAX = (MAX > r[v]) ? MAX: r[v];
		}
	}
        max(a, MAX, &MAXIMUM);
        if(a->num == 0) printf("Residual = %e %u %u %u Time = %e\n", MAXIMUM, n, m, p, time);
	delete [] r;
	delete [] id;
	delete [] id2;
	delete [] m3;
	delete [] m4;
}
void ERROR(arg *a, int *error)
{
    static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    int p = a->p;
    static int current_error = 0;
    static int t_in = 0;
    static int t_out = 0;
    pthread_mutex_lock(&mut);
    if(a->error != 0)
    {
        current_error = -1;
    }
    t_in++;
    if(t_in >= p)
    {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }
    else
    {
        while(t_in < p)
        {
            pthread_cond_wait(&c_in, &mut);
        }
    }
    *error = current_error;
    t_out++;
    if(t_out >= p)
    {
        current_error = 0;
        t_in = 0;
        pthread_cond_broadcast(&c_out);
    }
    else
    {
        while(t_out < p)
        {
            pthread_cond_wait(&c_out, &mut);
        }
    }
    pthread_mutex_unlock(&mut);
}
void* thread_func(void* args)
{
	arg* a = (arg*)args;
        int n = a->n;
        int m = a->m;
        int error = 0;
	double time = 0;
	double time2 = 0;
	double time_thread = 0;
	double norm = 0;
	int cpuNum  = get_nprocs();
        int* temp2;
        int* temp;
        double* block;
        double* help1;
        double* inv_block;
        block = new double[m*m];
        inv_block = new double[m*m];
        help1 = new double[m*m];
        temp = new int[m*m];
        temp2 = new int[n/m];
        for(int i = 0; i < n/m; i++)
        {
                temp2[i]=i;
        }
	cpu_set_t cpu;
	CPU_ZERO(&cpu);
	CPU_SET(cpuNum - a->num - 1, &cpu);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);
	matrixThreadsZeros(a);
	if(matrixThreadsData(a) != 0)
	{
		a->error = -2;
                delete [] block;
                delete [] inv_block;
                delete [] temp;
                delete [] help1;
                delete [] temp2;
                return 0;
	}
        matrixThreadsNorm(a, &norm);
        time = getTime();
        time_thread = getCPUTime();
        matrixThreadsInverse(a, norm, block, help1, temp, temp2, inv_block);
        a->time =  getCPUTime() - time_thread;
        time2 = getTime() - time;
        matrixThreadsData2(a);
        ERROR(a, &error);
        if(error == 0) discrepansyThread(a, time2);
        if(a->num == 0 && a->error == 0)
        {
                cout<<"General time = "<<time2<<"\n"<<endl;
        }
        delete [] block;
        delete [] inv_block;
        delete [] temp;
        delete [] help1;
        delete [] temp2;
        pthread_barrier_wait(a->barrier);
        return 0;
}
