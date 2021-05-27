all: ballAlg ballAlg-mpi

ballAlg: ballAlg.c 
	gcc -O3 -fopenmp ballAlg.c -o ballAlg -lm

ballAlg-mpi: ballAlg-mpi.c
	mpicc -O3 -fopenmp ballAlg-mpi.c -o ballAlg-mpi -lm
	
clean:
		rm -f *.o ballAlg ballAlg-mpi
