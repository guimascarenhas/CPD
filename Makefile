ballAlg: ballAlg.c 
	gcc -O3 -fopenmp ballAlg.c -o ballAlg -lm
	
clean:
		rm -f *.o ballAlg

valgrind:
		valgrind --leak-check=full ./ballAlg 2 5 0