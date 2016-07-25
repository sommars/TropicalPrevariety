#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"

using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

vector<vector<vector<GMP_Integer> > > CyclicN(int n) {
	vector<vector<vector<GMP_Integer> > > System;
	for (int i = 0; i != n-1; i++) {
		vector<vector<GMP_Integer> > Equation;
		for (int j = 0; j != n; j++) {
			vector<GMP_Integer> Monomial;
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

int main(int argc, char* argv[]) {
	if (argc != 2) {
		cout << "Internal error: expected an argument." << endl;
		return 1;
	}

	int n = atoi(argv[1]);

	vector<Hull> Hulls;
	vector<vector<vector<GMP_Integer> > > Cyc4 = CyclicN(n);
	vector<vector<vector<GMP_Integer> > >::iterator CycIt;
	for (CycIt=Cyc4.begin(); CycIt != Cyc4.end(); CycIt++) {
		Hulls.push_back(NewHull(*CycIt));
	}

	
	cout << "Hull count: " << Hulls.size() << endl;
	
	cout << "Found all hulls." << endl << endl << endl;
	
	
	vector<vector<C_Polyhedron> > Cones;
	vector<C_Polyhedron> Conevector;
	int HullIndex = 0;
	
	vector<Hull>::iterator it;
	for (it=Hulls.begin(); it != Hulls.end(); it++) {
		vector<Edge> Edges = (*it).Edges;
		vector<C_Polyhedron> HullCones;
		vector<Edge>::iterator itr;
		
		for (itr=Edges.begin(); itr != Edges.end(); itr++) {
			HullCones.push_back((*itr).Cone);
			//PrintCPolyhedron((*itr).Cone);cin.get();
		};
		
		if (HullIndex == 0) {
			Conevector = HullCones;
		} else {
			Cones.push_back(HullCones);
		};
		HullIndex++;
	};

	cout << "Conevector count: " << Conevector.size() << endl;
	cout << "Cones count: " << Cones.size() << endl;
	
	//Iterate through Cones
	vector<vector<C_Polyhedron> >::iterator ConesItr;
	int TreeLevel = 1;
	for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
		vector<C_Polyhedron> TestCones = *ConesItr;
		vector<C_Polyhedron> NewCones;
		//Iterate through Cones
		vector<C_Polyhedron>::iterator ConeItr;
		for (ConeItr=Conevector.begin(); ConeItr != Conevector.end(); ConeItr++) {
			//Iterate through TestCones
			C_Polyhedron Cone = *ConeItr;
			vector<C_Polyhedron>::iterator TestConesItr;
			for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
				C_Polyhedron NewCone = IntersectCones(*TestConesItr, Cone);
				NewCones.push_back(NewCone);
			};
		};
		Conevector = NewCones;
		printf("Finished level %d of tree with %d levels\n", TreeLevel, Cones.size());
		TreeLevel++;
	};
	
	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;
//	PrintCPolyhedrons(Conevector, false);
	
}
