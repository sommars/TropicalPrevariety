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
	
	cout << endl << endl << endl << endl << endl << endl << "Found all hulls." << endl << endl << endl;
	
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
		int IntersectionCount = 0;
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
					IntersectionCount++;
					if ( find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
						NewCones.push_back(NewCone);
					};
				};
			};
			ConeVector = NewCones;
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", TreeLevel, Cones.size(), ConeVector.size(),IntersectionCount);
			TreeLevel++;
		};
	}
	
	else if (string(argv[2]) == "new") {
		int ConeIntersectionCount = 0;
		int HullIndex = 0;	
		vector<Hull>::iterator it;
		for (it=Hulls.begin(); it != Hulls.end(); it++) {
			Hull H = (*it);
			vector<Edge> Edges = H.Edges;
			vector<Facet> Facets = H.Facets;
			map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap = H.PointToIndexMap;
			vector<vector<GMP_Integer> > Pts = H.Points;
			
			//Initialize ConeVector if it's the first convex hull.
			if (HullIndex == 0) {
				vector<Edge>::iterator itr;
				for (itr=Edges.begin(); itr != Edges.end(); itr++) {
					ConeVector.push_back((*itr).Cone);
				};
			} else {
				vector<C_Polyhedron> NewCones;
				vector<C_Polyhedron>::iterator ConeIterator;
				for (ConeIterator=ConeVector.begin();ConeIterator!=ConeVector.end();ConeIterator++) {
					//take a random vector from cone
					C_Polyhedron NewCone = *ConeIterator;

					vector<vector<GMP_Integer> > Rays = GeneratorSystemToPoints(NewCone.minimized_generators());
					//WARNING: this only works for cyclic n!
					GMP_Integer RandomArray[n];
					for (int i = 0; i < n; i++) {
						RandomArray[i] = 0;
					};
					vector<vector<GMP_Integer> >::iterator RayIterator;
					for (RayIterator=Rays.begin();RayIterator != Rays.end(); RayIterator++) {
						vector<GMP_Integer> Ray = *RayIterator;
						vector<GMP_Integer>::iterator RayValue;
						int RayIndex = 0;
						for (RayValue=Ray.begin();RayValue != Ray.end();RayValue++){
							RandomArray[RayIndex]+=*RayValue;
							RayIndex++;
						};
					};
					
					//WARNING: this only works for cyclic n!
					vector<GMP_Integer> RandomVector(RandomArray, + RandomArray + n);
					
					//take initial form using that vector.
					vector<vector<GMP_Integer> > InitialForm = FindInitialForm(Pts, RandomVector);
					
					set<GMP_Integer> InitialIndices;
					vector<vector<GMP_Integer> >::iterator InitialFormItrr;
					for (InitialFormItrr=InitialForm.begin(); InitialFormItrr != InitialForm.end(); InitialFormItrr++) {
						InitialIndices.insert(PointToIndexMap[*InitialFormItrr]);
					};
					//explore edge skeleton
					vector<GMP_Integer> EdgesToTest; //list of indices of edges to visit.
					vector<Edge>::iterator EdgeItr;
					GMP_Integer EdgeIndex = 0;
					for(EdgeItr = Edges.begin(); EdgeItr != Edges.end(); EdgeItr ++) {
						if (IntersectSets(InitialIndices, (*EdgeItr).PointIndices).size() > 0) {
							EdgesToTest.push_back(EdgeIndex);
						};
						EdgeIndex++;
					};
					
					vector<GMP_Integer> PretropGraphEdges;
					vector<GMP_Integer> NotPretropGraphEdges;
					
					while(!EdgesToTest.empty()) {
						GMP_Integer EToT = EdgesToTest.back();
						EdgesToTest.pop_back();
						GMP_Integer EInd = 0;
						vector<Edge>::iterator EdgeItr2;
						Edge EdgeToTest;
						for(EdgeItr2 = Edges.begin(); EdgeItr2 != Edges.end(); EdgeItr2++) {
							if (EInd == EToT) {
								EdgeToTest = *EdgeItr2;
								break;
							};
							EInd++;
						};
						C_Polyhedron TempCone = IntersectCones(EdgeToTest.Cone, NewCone);
						ConeIntersectionCount++;
						//Note: this is why we want to use reduced cyclic n. Consider switching this to if it has any rays.
						if (TempCone.affine_dimension() > 0) {
							PretropGraphEdges.push_back(EToT);
							
							if ( find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
								NewCones.push_back(NewCone);
							};	
							set<GMP_Integer>::iterator NeighborItr;
							for(NeighborItr=EdgeToTest.NeighborIndices.begin();NeighborItr!=EdgeToTest.NeighborIndices.end(); NeighborItr++){
								GMP_Integer Neighbor = *NeighborItr;
								//if Neighbor not in PretropGraphEdges:
								if ( find(PretropGraphEdges.begin(), PretropGraphEdges.end(), Neighbor) == PretropGraphEdges.end() ) {
									if ( find(NotPretropGraphEdges.begin(), NotPretropGraphEdges.end(), Neighbor) == NotPretropGraphEdges.end() ) {
										if ( find(EdgesToTest.begin(), EdgesToTest.end(), Neighbor) == EdgesToTest.end() ) {
											EdgesToTest.push_back(Neighbor);
										};
									};
								};
								EdgeToTest.NeighborIndices;
							};
						} else {
							NotPretropGraphEdges.push_back(EToT);
						};
					}

					
					
				};
				ConeVector = NewCones;
			};
			HullIndex++;
		};

		cout << "ConeVector count: " << ConeVector.size() << endl;

	
	
	
	};
	
	
	
	
	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;
	vector<Generator> gv;
	vector<C_Polyhedron>::iterator PolyItr;
	vector<vector<vector<GMP_Integer> > > InitialForms;
	for (PolyItr=ConeVector.begin(); PolyItr != ConeVector.end(); PolyItr++) {
		Generator_System gs = (*PolyItr).minimized_generators();
		for (Generator_System::const_iterator gsi = gs.begin(),
		gs_end = gs.end(); gsi != gs_end; ++gsi) {
			Generator gen = *gsi;
			if (gen.is_point() or gen.is_line()) {
				continue;
			}
			else if ( find(gv.begin(), gv.end(), gen) == gv.end() ) {
				
				vector<GMP_Integer> Pt = GeneratorToPoint(gen);
				vector<Hull>::iterator HullItr;
				vector<vector<GMP_Integer> > MyInitialForms;
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
	//				PrintPoints(InForm);
				//	cout << endl;
					sort(PtIndices.begin(),PtIndices.end());
					MyInitialForms.push_back(PtIndices);
				}
				if (find(InitialForms.begin(), InitialForms.end(), MyInitialForms) == InitialForms.end() ) {
					InitialForms.push_back(MyInitialForms);
					gv.push_back(gen);
					cout << "Pretropism: " << gen << endl;
				}
			};
		};
	};
	cout << "Number of pretropisms found: " << gv.size() << endl;
}
