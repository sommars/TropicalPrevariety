#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"
#include <algorithm>

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
	if (argc != 3) {
		cout << "Internal error: expected two arguments." << endl;
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
	
	vector<C_Polyhedron> ConeVector;	
	
	if (string(argv[2]) == "refinement") {
		vector<vector<C_Polyhedron> > Cones;

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
				ConeVector = HullCones;
			} else {
				Cones.push_back(HullCones);
			};
			HullIndex++;
		};

		cout << "ConeVector count: " << ConeVector.size() << endl;
		cout << "Cones count: " << Cones.size() << endl;
	
		//Iterate through Cones
		vector<vector<C_Polyhedron> >::iterator ConesItr;
		int TreeLevel = 1;
		for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
			vector<C_Polyhedron> TestCones = *ConesItr;
			vector<C_Polyhedron> NewCones;
			//Iterate through Cones
			vector<C_Polyhedron>::iterator ConeItr;
			for (ConeItr=ConeVector.begin(); ConeItr != ConeVector.end(); ConeItr++) {
				//Iterate through TestCones
				C_Polyhedron Cone = *ConeItr;
				vector<C_Polyhedron>::iterator TestConesItr;
				for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
					C_Polyhedron NewCone = IntersectCones(*TestConesItr, Cone);
					//PrintPolyhedron(NewCone);
					if ( find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
						NewCones.push_back(NewCone);
					};
				};
			};
			ConeVector = NewCones;
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level\n", TreeLevel, Cones.size(), ConeVector.size());
			TreeLevel++;
		};
	}
	else if (string(argv[2]) == "new") {
	
	
	
	};
	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;	
	cout << argv[2] << endl;
	vector<Generator> gv;
	vector<C_Polyhedron>::iterator PolyItr;
	for (PolyItr=ConeVector.begin(); PolyItr != ConeVector.end(); PolyItr++) {
		Generator_System gs = (*PolyItr).minimized_generators();
		for (Generator_System::const_iterator gsi = gs.begin(),
		gs_end = gs.end(); gsi != gs_end; ++gsi) {
			Generator gen = *gsi;
			if (gen.is_point()) {
				continue;
			}
			else if ( find(gv.begin(), gv.end(), gen) == gv.end() ) {
				gv.push_back(gen);
				cout << "Pretropism: " << gen << endl;
				
				vector<GMP_Integer> Pt = GeneratorToPoint(gen);
				vector<Hull>::iterator HullItr;
				for (HullItr=Hulls.begin(); HullItr != Hulls.end(); HullItr++) {
					Hull MyHull = (*HullItr);
					vector<vector<GMP_Integer> > TestPts = MyHull.Points;
					vector<vector<GMP_Integer> > InForm = FindInitialForm(TestPts, Pt);
					if (InForm.size() < 2) {
						cout << "INTERNAL ERROR: Incorrect pretropism found" << endl;
						PrintPoints(InForm);
						cin.get();
					};

					vector<GMP_Integer> PtIndices;
			//		cout << InForm.size() << endl;
					
	
					vector<vector<GMP_Integer> >::iterator InFormIter;
					for (InFormIter=InForm.begin(); InFormIter!=InForm.end(); InFormIter++) {
						PtIndices.push_back(MyHull.PointToIndexMap[*InFormIter]);
					}
		//			PrintPoint(PtIndices);
	//				PrintPoints(InForm);
				//	cout << endl;
				}
			};
		};
	};
	cout << "Number of pretropisms found: " << gv.size() << endl;
}
