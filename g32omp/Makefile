all: ballAlg ballAlg-omp

ballAlg: ballAlg.c 
	gcc -O3 -fopenmp ballAlg.c -o ballAlg -lm

ballAlg-omp: ballAlg-omp.c
	gcc -O3 -fopenmp ballAlg-omp.c -o ballAlg-omp -lm
	
clean:
		rm -f *.o ballAlg ballAlg-omp
