
//   echo 'localhost slots=25' > hostfile
//   mpicc barrier_fin.c -lstdc++ -o barrier_fin
//   mpirun -c 8 --hostfile hostfile barrier_fin



#include <stdio.h>
#include <mpi.h>


void barrier(int rank, int size)
{
    int power = 1;
    int buf = 0;
    
    
    // 1
    while (power < size)
    {
        power *= 2;
        if (rank % power)
        {
            int RecvRank = rank - power / 2;
            if (RecvRank >= 0)
            {
                MPI_Send(&buf, 1, MPI_INT, RecvRank, 0, MPI_COMM_WORLD);
                break;
            }
        }
        else
        {
            int SendRank = rank + power / 2;
            if (SendRank  < size)
            {
                MPI_Recv(&buf, 1, MPI_INT,  SendRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("1       %d   =>   %d\n", SendRank, rank);
            }
                
                
        }
       
    }
    
    
    // 2
    while (power > 1)
    {
        if (rank % power)
        {
            int SendRank = rank - power / 2;
            if (SendRank >= 0)
            {
                MPI_Recv(&buf, 1, MPI_INT, SendRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printf("2       %d   =>   %d\n", SendRank, rank);
            }
        }
        else
        {
            int RecvRank = rank + power / 2;
            if (RecvRank < size)
            {
                MPI_Send(&buf, 1, MPI_INT, RecvRank, 1, MPI_COMM_WORLD);
            }
        }
        
        power /= 2;
    }

}




int main(int argc, char *argv[])
{
    int ProcNum, ProcRank, RecvRank;
    MPI_Status Status;
    double time1, time2;
    
    MPI_Initialized(&argc);
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    
    
    time1 = MPI_Wtime();
    
    if(ProcNum <= 1)
    {
        printf("Please define some proceses");
        MPI_Finalize();
        return 0;
    }
    
    
    
    barrier(ProcRank, ProcNum );
    
    
    time2 = MPI_Wtime() - time1;
    
    //printf("time = %f \n", time2);
    
    MPI_Finalize();
    return 0;
}


