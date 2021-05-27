/*
 *   Send / Recv
 *
 *   Jos� Monteiro, DEI / IST
 *
 *   Last modification: 2 November 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main (int argc, char *argv[]) {

    MPI_Status status;
    int id, p,
	i, rounds;
    double secs;
    char message[20] = "Random message";

    MPI_Init (&argc, &argv);

    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &p);

    rounds = atoi(argv[1]);

    MPI_Barrier (MPI_COMM_WORLD);
    secs = - MPI_Wtime();

    for(i = 0; i < rounds; i++){
	if(!id){
	    MPI_Send(&i, 1, MPI_INT, 1, i, MPI_COMM_WORLD);
	    MPI_Recv(&i, 1, MPI_INT, p-1, i, MPI_COMM_WORLD, &status);
        printf("id: %d, p: %d  from rank %d, with tag %d and error code %d.", id,p,status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);
	}
	else{
	    MPI_Recv(&i, 1, MPI_INT, id-1, i, MPI_COMM_WORLD, &status);
	    MPI_Send(&i, 1, MPI_INT, (id+1)%p, i, MPI_COMM_WORLD);
        printf("id: %d, p: %d  from rank %d, with tag %d and error code %d.", id,p,status.MPI_SOURCE, status.MPI_TAG, status.MPI_ERROR);
	}
    }

    MPI_Barrier (MPI_COMM_WORLD);
    secs += MPI_Wtime();

    if(!id){
	printf("Rounds= %d, N Processes = %d, Time = %12.6f sec,\n",
	       rounds, p, secs);
	printf ("Average time per Send/Recv = %6.2f us\n",
		secs * 1e6 / (2*rounds*p));
   }
   MPI_Finalize();
   return 0;
}
