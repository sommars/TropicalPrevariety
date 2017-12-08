main:
	clear
	/bin/rm -f *.o *.out
	g++ -O3 -std=c++11 -c prevariety_util.cpp
	g++ -O3 -std=c++11 -c cone_tree.cpp
	g++ -O3 -std=c++11 -g cone_tree.o prevariety_util.o cone_intersection.cpp -lppl -lgmpxx -lgmp -o prevariety.out
