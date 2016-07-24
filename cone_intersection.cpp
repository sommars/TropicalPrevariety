#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"

using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

list<list<list<GMP_Integer> > > CyclicN(int n) {
	list<list<list<GMP_Integer> > > System;
	for (int i = 0; i != n-1; i++) {
		list<list<GMP_Integer> > Equation;
		for (int j = 0; j != n; j++) {
			list<GMP_Integer> Monomial;
			set<GMP_Integer> OneSet;
			for (int k = 0; k != i+1; k++) {
				OneSet.insert((j+k)%n);
			};
			for (int l = 0; l != n; l++) {
				if (OneSet.find(l) != OneSet.end()) {
					Monomial.push_back(1);
				} else {
					Monomial.push_back(0);
				};
			};
			Equation.push_back(Monomial);
		};
		System.push_back(Equation);
	};
	return System;
}

int main() {
	list<Hull> Hulls;
	list<list<list<GMP_Integer> > > Cyc4 = CyclicN(4);
	list<list<list<GMP_Integer> > >::iterator CycIt;
	for (CycIt=Cyc4.begin(); CycIt != Cyc4.end(); CycIt++) {
		Hulls.push_back(NewHull(*CycIt));
	}

	
	cout << "Hull count: " << Hulls.size() << endl;
	
	cout << "Found all hulls." << endl << endl << endl;
	
	
	list<list<C_Polyhedron> > Cones;
	list<C_Polyhedron> ConeList;
	int HullIndex = 0;
	
	list<Hull>::iterator it;
	for (it=Hulls.begin(); it != Hulls.end(); it++) {
		list<Edge> Edges = (*it).Edges;
		list<C_Polyhedron> HullCones;
		list<Edge>::iterator itr;
		for (itr=Edges.begin(); itr != Edges.end(); itr++) {
			PrintCPolyhedron((*itr).Cone);
			HullCones.push_back((*itr).Cone);
		};
		
		if (HullIndex == 0) {
			ConeList = HullCones;
		} else {
			Cones.push_back(HullCones);
		};
		HullIndex++;
	};
	
	cout << "ConeList count: " << ConeList.size() << endl;
	cout << "Cones count: " << Cones.size() << endl;
	
	//Iterate through Cones
	list<list<C_Polyhedron> >::iterator ConesItr;
	for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
		list<C_Polyhedron> TestCones = *ConesItr;
		list<C_Polyhedron> NewCones;
		//Iterate through Cones
		list<C_Polyhedron>::iterator ConeItr;
		for (ConeItr=ConeList.begin(); ConeItr != ConeList.end(); ConeItr++) {
			//Iterate through TestCones
			C_Polyhedron Cone = *ConeItr;
			list<C_Polyhedron>::iterator TestConesItr;
			for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
				C_Polyhedron NewCone = IntersectCones(*TestConesItr, Cone);
				NewCones.push_back(NewCone);
			};
		};
		ConeList = NewCones;
	};
	
	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;
//	PrintCPolyhedrons(ConeList, false);
	
}
