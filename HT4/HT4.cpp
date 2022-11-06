#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <time.h> 

int main(int argc, char** argv) {
    int ProcRank, ProcNum, N = 1000000, CountIn = 0;
    MPI_Status Status;
    double* x = new double[N];
    double* y = new double[N];
    srand(time(NULL));
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double* RecievedPi = new double[ProcNum];
    double* RecievedSec = new double[ProcNum];
    int* RecievedNum = new int[ProcNum];
    double* RecievedStatus = new double[ProcNum];
    if (ProcRank == 0){
        for (int i = 0; i < N; i++){
            x[i] = (double)rand() / RAND_MAX;
            y[i] = (double)rand() / RAND_MAX; 
        }
    }
    MPI_Bcast(x, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int k;
    if (ProcNum == 1){
        k = N / (ProcNum);
    }
    else {
        k = N / (ProcNum - 1);
    }
    int i2 = k * ProcRank;
    if (ProcRank == ProcNum - 1) i2 = N;

    if (ProcRank != 0){
        double start = MPI_Wtime();
        for (int i = 0; i < i2; i++){
            double z = sqrt(x[i] * x[i] + y[i] * y[i]);
            if (z <= 1){
                CountIn++;
            }
        }
        double pi = 4 * ((double)CountIn / (double)i2);
        double end = MPI_Wtime();
        double seconds = (end - start);
        MPI_Send(&pi, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&seconds, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&i2, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
    }
    else if (ProcRank == 0){
        for (int i = 1; i < ProcNum; i++){
            MPI_Recv(&RecievedPi[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&RecievedSec[i], 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&RecievedNum[i], 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &Status);
            RecievedStatus[i] = Status.MPI_SOURCE;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (ProcRank == 0) {
        for (int i = 1; i < ProcNum; i++) {
            std::cout << "process: " << RecievedStatus[i] << ", pi value: " << RecievedPi[i] << ", seconds: " << RecievedSec[i]
                << ", number of dots: " << RecievedNum[i] << std::endl;
        }
    }
    MPI_Finalize();

    delete[] x;
    delete[] y;
    delete[] RecievedPi;
    delete[] RecievedSec;
    delete[] RecievedNum;
    delete[] RecievedStatus;
    return 0;
}
