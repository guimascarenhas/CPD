ballAlg: ballAlg.c 
	gcc -O3 -fopenmp ballAlg.c -o ballAlg -lm
	
clean:
		rm -f *.o ballAlg