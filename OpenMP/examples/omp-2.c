/*
DESCRIPTION:
		Parallelizing an inner loop with dependences

		for (iter=0; iter<numiter; iter++) {
			for (i=0; i<size-1; i++) {
				V[i] = f( V[i], V[i+1] );
			}
		}

**************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

#define TOTALSIZE 100
#define NUMITER 2000

/*
* DUMMY FUNCTION
*/
#define f(x,y)	((x+y)/2.0)


/* MAIN: PROCESS PARAMETERS */
int main(int argc, char *argv[]) {

  /* VARIABLES */
  int i, iter;

  /* DECLARE VECTOR AND AUX DATA STRUCTURES */
  double *V = (double *) malloc(TOTALSIZE * sizeof(double));
  double *aux = (double *) malloc(TOTALSIZE * sizeof(double));

  /* 1. INITIALIZE VECTOR */
  for(i = 0; i < TOTALSIZE; i++) {
    V[i]= 0.0 + i;
  }

    double exec_time;
    exec_time = -omp_get_wtime();
  /* 2. ITERATIONS LOOP */
  for(iter = 0; iter < NUMITER; iter++) {
    if(V[0] != (NUMITER-1)){
    #pragma omp parallel private(i)
    {
    #pragma omp for 
    for(i = 0; i < TOTALSIZE; i++) {
     aux[i] = V[i];
    }
    /* 2.2. END ITERATIONS LOOP */



    /* 2.1. PROCESS ELEMENTS */

    #pragma omp for 
    for(i = 0; i < TOTALSIZE-1; i++) {
      V[i] = f(aux[i], aux[i+1]);
    }


    /* 2.2. END ITERATIONS LOOP */
    }
    }
  }
  
  exec_time += omp_get_wtime();
  /* 3. OUTPUT FINAL VALUES */
  printf("Output:\n"); 
  for(i = 0; i < TOTALSIZE; i++) {
    printf("%4d %f\n", i, V[i]);
  }
  printf("Exec time: %f\n", exec_time);
}
