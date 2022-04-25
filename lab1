#include <iostream>
#include "mpi.h"


int main(int argc, char *argv[])
{
    int ProcNum, ProcRank, M, RecvRank;
    MPI_Status Status;
    M = 10;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    for (int i = 0; i < M; i++) {
        if (ProcRank == 0 && i == 0) {
            MPI_Send(&ProcRank, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
            continue;
        }

        if (ProcRank == 0)
            MPI_Recv(&RecvRank, 1, MPI_INT, ProcNum - 1, 1, MPI_COMM_WORLD, &Status);
        else
            MPI_Recv(&RecvRank, 1, MPI_INT, ProcRank-1, 1, MPI_COMM_WORLD, &Status);

        std::cout << ProcRank << std::endl;

        if(ProcRank == ProcNum-1)
            MPI_Send(&ProcRank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        else
            MPI_Send(&ProcRank, 1, MPI_INT, ProcRank+1, 1, MPI_COMM_WORLD);
    }
    MPI_Finalize();

    return 0;
}
