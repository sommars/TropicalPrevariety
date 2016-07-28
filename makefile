main:
	clear
	/bin/rm -f *.o *.out
	g++ -O3 -c prevariety_util.cpp
	g++ -O3 -g prevariety_util.o cone_intersection.cpp -lppl -lgmpxx -lgmp -o prevariety.out
