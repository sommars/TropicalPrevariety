#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>

using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
int InnerProduct(list<int> V1, list<int> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
	};
	int Result = 0;
	list<int>::iterator it1;
	list<int>::iterator it2;
	it2 = V2.begin();
	for (it1=V1.begin(); it1 != V1.end(); it1++) {
		Result = Result + (*it1) * (*it2);
		it2++;
	}
	return Result;
}


//------------------------------------------------------------------------------
list<list<int> > FindInitialForm(list<list<int> > Points, list<int> Vector) {
	/* 
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};

	list<list<int> > IF;
	list<list<int> >::iterator itr;

	itr=Points.begin();
	IF.push_back(*itr);
	int MinimalIP = InnerProduct(Vector, *itr);
	itr++;

	for (itr; itr != Points.end(); itr++) {
		int IP = InnerProduct(Vector, *itr);
		if (MinimalIP > IP) {
			MinimalIP = IP;
			IF.clear();
			IF.push_back(*itr);
		} else if (IP == MinimalIP) {
			IF.push_back(*itr);
		}
	}
	
	return IF;
}

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(list<list<int> > Points) {
	Generator_System gs;
	list<list<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		list<int>Point=*itr;
		list<int>::iterator it;
		Linear_Expression LE;
		int VarIndex = 0;
		for (it=Point.begin(); it != Point.end(); it++) {
			LE = LE + Variable(VarIndex) * (*it);
			VarIndex++;
		}
		gs.insert(point(LE));
	}
	C_Polyhedron ph = C_Polyhedron(gs);
	return ph;
}
