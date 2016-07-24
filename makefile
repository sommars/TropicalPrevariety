main:
	clear
	/bin/rm -f *.o *.out
	g++ -c prevariety_util.cpp
	g++ -g prevariety_util.o cone_intersection.cpp -lppl -lgmpxx -lgmp -o prevariety.out
