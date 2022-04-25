#include <iostream>
#include "mpi.h"


int main(int argc, char *argv[])
{
    int ProcNum, ProcRank, M, RecvRank;
    int SendBuf;
    MPI_Status Status;
    M = 5;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    SendBuf = 0;
    MPI_Bcast(&SendBuf, 1, MPI_INT, SendBuf, MPI_COMM_WORLD);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < ProcNum; j++)
        {
            if(SendBuf == ProcRank-1 || SendBuf==ProcNum-1 && ProcRank==0)
                std::cout << SendBuf << std::endl;
            SendBuf++;
            if (SendBuf == ProcNum)
                SendBuf = 0;
            MPI_Bcast(&SendBuf, 1, MPI_INT, SendBuf, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();

    return 0;
}
