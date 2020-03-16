#include"Jordan_Header.hpp"
using namespace std;
int main(int argc, char* argv[])
    /*  argv[1] = matrix_size -- размер матрицы;
        argv[2] = block_size -- размер блока;
        argv[3] = number_of_threads -- число потоков;
        Значение argc должно быть не меньше 4: a.out + n + m + p; */
{
    if (argc < 4)
    {
            cout<<"Usage: "<<argv[0]<<" <n> <m> <p> <file_name>"<<endl;
            return -1;
    }
    if (atoi(argv[1]) < atoi(argv[2]))
    {
            cout<<"Введенный размер блока больше размера матрицы."<<endl;
            return -1;
    }
    if (atoi(argv[3]) < 1)
    {
            cout<<"Некорректное количество потоков"<<endl;
            return -1;
    }
    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int p = atoi(argv[3]);
    double *matrix1, *matrix2;
    matrix1 = new double[n*n];
    matrix2 = new double[n*n];
    //row = new double[m];
    if(!matrix1)
    {
        cout<<"Not enough memory!"<<endl;
        return -1;
    }
    if(!matrix2)
    {
        cout<<"Not enough memory!"<<endl;
        delete [] matrix1;
        return -1;
    }
    pthread_barrier_t barrier;
    pthread_t  *threads;
    threads = new pthread_t[p];
    if(!threads)
    {
        cout<<"Cannot create array of threads!"<<endl;
        return 2;
    }
    arg* a;
    a = new arg[p];
	if(!a)
    {
        cout<<"Not enough memory!"<<endl;
        delete[] threads;
        return 2;
    }
    pthread_barrier_init(&barrier, 0, p);
    for (int k = 0; k < p; k++)
    {
            a[k].num = k;
            a[k].p = p;
            a[k].n = n;
            a[k].m = m;
            a[k].a = matrix1;
            a[k].b = matrix2;
            a[k].barrier = &barrier;
            a[k].error = 0;
            if(argc == 5)
            {
                    a[k].file_name = argv[4];
            }
            else
            {
                    a[k].file_name = NULL;
            }
    }
    for(int k = 1; k < p; k++)
    {
        if(pthread_create(&threads[k], NULL, thread_func, a + k) != 0)
        {
            cout<<"Cannot create thread"<<k<<endl;
            pthread_barrier_destroy(&barrier);
            delete[] matrix1;
            delete[] matrix2;
            delete[] threads;
            delete[] a;
            return -1;
        }
    }
    thread_func(a + 0);
    for(int k = 0; k < p; k++)
    {
        if(a[k].error != 0) 
        {
            cout<<"Error!"<<endl;
            pthread_barrier_destroy(&barrier);
            delete[] matrix1;
            delete[] matrix2;
            delete[] threads;
            delete[] a;
            return -1;
        }
    }
    for(int k = 0; k < p; k++)
    {
        cout<<"Time of thread number "<<a[k].num<<" = "<<a[k].time<<endl;
    }
    cout<<endl;
    printMatrix(matrix2, n, m, (n > 10)? 10: n);
    delete[] matrix1;
    delete[] matrix2;
    pthread_barrier_destroy(&barrier);
    delete[] threads;
    delete[] a;
    //delete[] row;
    exit(0);
}
