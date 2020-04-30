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
    double norm = 0, residual;
    int n, m, whole, max_rows;
    char* file_name = NULL;
    if(argc == 3 || argc == 4)
    {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        if(n < m)
        {
            if(my_rank == 0) cout<<"\nIncorrect size of block!\n"<<endl;
            MPI_Finalize();
            return 0;
        }
    }
    else
    {
        if(my_rank == 0) cout<<"\nUsage: "<<argv[0]<<" <n> <m> <file_name>\n"<<endl;
        MPI_Finalize();
        return 0;
    }
    whole = n/m;
    max_rows = whole/p + 1;
    matrix = new double[n*m*max_rows];
    id = new double[n*m*max_rows];
    row_a = new double[n*m];
    row_b = new double[n*m];
    if(matrix == NULL)
    {
        cout<<"\nNot enough memory!\n"<<endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(id == NULL)
    {
        cout<<"\nNot enough memory!\n"<<endl;
        delete [] matrix;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(row_a == NULL)
    {
        printf("\nNot enough memory!\n");
        delete [] matrix;
        delete [] id;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(row_b == NULL)
    {
        printf("\nNot enough memory!\n");
        delete [] matrix;
        delete [] id;
        delete [] row_a;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }
    if(createMatrix(matrix, row_a, n, m, p, my_rank, argc, file_name) == -1)
    {
        delete [] matrix;
        delete [] id;
        delete [] row_a;
        delete [] row_b;
        MPI_Finalize();
        return -1;
    }
    createMatrixFormula(id, n, m, p, my_rank, 0);
    printMatrixMPI(matrix, my_rank, p, n, m, row_a);
    // printMatrixMPI(id, my_rank, p, n, m, row_a);
    MPINorm(matrix, row_a, n, m, my_rank, p, &norm);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_of_start = 0.0, total_time = 0.0;
    if(my_rank == 0) time_of_start = MPI_Wtime();
    if(matrixInverseMPI(matrix, id, row_a, row_b, norm, n, m, p, my_rank) == -1)
    {
        if(my_rank == 0) cout<<"\nМетод неприменим для данного m!\n"<<endl;
        delete [] matrix;
        delete [] id;
        delete [] row_a;
        delete [] row_b;
        MPI_Abort(MPI_COMM_WORLD, -2);
        return -1;
    }
    if(my_rank == 0) total_time = MPI_Wtime() - time_of_start;
    MPI_Barrier(MPI_COMM_WORLD);
    printMatrixMPI(id, my_rank, p, n, m, row_a);
//    makeRowsMatrix(matrix, id, n, m, my_rank, p);
//    if(argc == 4)
//    {
//        file_name = argv[3];
//        if(readMatrixFromFile(id, row_a, n, m, p, my_rank, file_name) == -1)
//        {
//            delete [] matrix;
//            delete [] id;
//            delete [] row_a;
//            delete [] row_b;
//            MPI_Finalize();
//            return -1;
//        }
//    }
//    else
//    {
//        createMatrix(id, n, m, p, my_rank, TYPE_OF_MATRIX);
//    }
    residual = residualMPI(matrix, id, row_a, n, m, my_rank, p, max_rows, argc, file_name);
    if(my_rank == 0) cout<<"Residual = "<<setw(SETW_CONSTANT)<<residual<<", time = "<<setw(SETW_CONSTANT)<<total_time<<", n = "<<setw(5)<<n<<", m = "<<setw(5)<<m<<", p = "<<setw(5)<<p<<endl;
    delete [] matrix;
    delete [] id;
    delete [] row_a;
    delete [] row_b;
    MPI_Finalize();
    return 0;
}
