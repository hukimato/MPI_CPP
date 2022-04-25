#include <iostream>
#include "mpi.h"


int main(int argc, char *argv[])
{
    int ProcNum, ProcRank;
    int M = 5;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int* SendBuf = new int[ProcNum];
    for (int i = 0; i < ProcNum; i++)
    {
        SendBuf[i] = ProcRank;
    }
    int RecvBuf;

    MPI_Scatter(SendBuf, 1, MPI_INT,
        &RecvBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < ProcNum; j++)
        {
            if(RecvBuf == ProcRank-1 || RecvBuf == ProcNum-1 && ProcRank == 0)
                std::cout << RecvBuf << std::endl;
            
            RecvBuf++;
            if (RecvBuf == ProcNum)
                RecvBuf = 0;

            MPI_Scatter(SendBuf, 1, MPI_INT,
                &RecvBuf, 1, MPI_INT, RecvBuf, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();

    return 0;
}
