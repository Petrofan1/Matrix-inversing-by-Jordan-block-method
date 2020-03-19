#include "Zhordan_Header_MPI.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    if(MPI_Init(&argc, &argv) != 0)
    {
        return -1;
    }
    int my_rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    double* matrix;
    double* id;
    double* row_a;
    double* row_b;
    double norm, residual;
    int n, m, whole, max_rows;
    char* file_name;
    if(argc == 3 || argc == 4)
    {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        if(n < m)
        {
            cout<<"Incorrect size of block!"<<endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return -1;
        }
    }
    else
    {
        cout<<"Usage: "<<endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    whole = n/m;
    max_rows = whole/p + 1;
    matrix = new double[n*m*max_rows];
    id = new double[n*m*max_rows];
    row_a = new double[n*m];
    row_b = new double[n*m];
    if(matrix == NULL)
    {
        cout<<"Not enough memory!"<<endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(id == NULL)
    {
        cout<<"Not enough memory!"<<endl;
        delete [] matrix;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(row_a == NULL)
    {
        printf("Not enough memory!");
        delete [] matrix;
        delete [] id;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(row_b == NULL)
    {
        printf("Not enough memory!");
        delete [] matrix;
        delete [] id;
        delete [] row_a;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(argc == 4)
    {
        file_name = argv[3];
        if(readMatrixFromFile(matrix, row_a, n, m, p, my_rank, file_name) == -1)
        {
            delete [] matrix;
            delete [] id;
            delete [] row_a;
            delete [] row_b;
            MPI_Finalize();
            return -1;
        }
    }
    else
    {
        createMatrix(matrix, n, m, p, my_rank, TYPE_OF_MATRIX);
    }
    createMatrix(id, n, m, p, my_rank, 0);
    printMatrixMPI(matrix, my_rank, p, n, m, row_a);
    MPINorm(matrix, row_a, n, m, my_rank, p, &norm);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_of_start = 0.0, total_time = 0.0;
    if(my_rank == 0)
    {
        time_of_start = MPI_Wtime();
    }
    if(matrixInverseMPI(matrix, id, row_a, row_b, norm, n, m, p, my_rank) == -1)
    {
        if(my_rank == 0)
        {
            cout<<"Метод неприменим для данного m!"<<endl;
        }
        delete [] matrix;
        delete [] id;
        delete [] row_a;
        delete [] row_b;
        MPI_Abort(MPI_COMM_WORLD, -2);
        return -1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0)
    {
        total_time = MPI_Wtime() - time_of_start;
    }
    printMatrixMPI(id, my_rank, p, n, m, row_a);
    if(my_rank == 0)
    {
        cout<<"Time = "<<total_time<<endl;
    }
    residual = residualMPI();
    cout<<"Residual = "<<residual<<endl;
    delete [] matrix;
    delete [] id;
    delete [] row_a;
    delete [] row_b;
    MPI_Finalize();
    return 0;
}
