#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"
#include <algorithm>
#include "tbb/tbb.h"
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

double ListAppendTime;
double ContainmentTime;
int ConeIntersectionCount;
vector<vector<vector<GMP_Integer> > > CyclicN(int n, bool Reduced) {
	vector<vector<vector<GMP_Integer> > > System;
	int Length = n;
	if (Reduced == true) {
		Length = n - 1;
	};
	
	for (int i = 0; i != n-1; i++) {
		vector<vector<GMP_Integer> > Equation;
		for (int j = 0; j != n; j++) {
			vector<GMP_Integer> Monomial;
			set<GMP_Integer> OneSet;
			for (int k = 0; k != i+1; k++) {
				OneSet.insert((j+k)%n);
			};
			for (int l = 0; l != Length; l++) {
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

//------------------------------------------------------------------------------
vector<C_Polyhedron> TraverseEdgeSkeleton(C_Polyhedron NewCone, vector<vector<GMP_Integer> > Pts, map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap, vector<Edge> Edges) {

	vector<C_Polyhedron> Result;
	
	//take a random vector from cone
	vector<vector<GMP_Integer> > Rays = GeneratorSystemToPoints(NewCone.minimized_generators());
	vector<GMP_Integer> RandomVector(Rays[0].size(),0);
	for (int i = 0; i != Rays.size(); i++) {
		vector<GMP_Integer> Ray = Rays[i];
		for (int j = 0; j != Ray.size(); j++) {
			RandomVector[j] += Ray[j];
		};
	};
	
	//take initial form using that vector.
	vector<vector<GMP_Integer> > InitialForm = FindInitialForm(Pts, RandomVector);
	
	set<GMP_Integer> InitialIndices;
	for (vector<vector<GMP_Integer> >::iterator InitialFormItr=InitialForm.begin();
		InitialFormItr != InitialForm.end(); InitialFormItr++)
	{
		InitialIndices.insert(PointToIndexMap[*InitialFormItr]);
	};
	
	// List of indices of edges to visit.
	vector<int> EdgesToTest; 
	for(int EdgeIndex = 0; EdgeIndex != Edges.size(); EdgeIndex++) {
		if (IntersectSets(InitialIndices, Edges[EdgeIndex].PointIndices).size() > 0) {
			EdgesToTest.push_back(EdgeIndex);
		};
	};
	
	// Explore edge skeleton
	set<int> PretropGraphEdges;
	set<int> NotPretropGraphEdges;
	while(!EdgesToTest.empty()) {
		int EdgeToTestIndex = EdgesToTest.back();
		EdgesToTest.pop_back();
		Edge EdgeToTest = Edges[EdgeToTestIndex];
		clock_t ContainmentBegin = clock();
		bool ConeDoesContain = EdgeToTest.Cone.contains(NewCone);
		ContainmentTime += double(clock() - ContainmentBegin);

		if (ConeDoesContain) {
			Result.push_back(NewCone);
			return Result;
		};
		C_Polyhedron TempCone = IntersectCones(EdgeToTest.Cone, NewCone);
		ConeIntersectionCount++;

		if (TempCone.affine_dimension() > 0) {
			PretropGraphEdges.insert(EdgeToTestIndex);
			
			clock_t begin = clock();
			if ( find(Result.begin(), Result.end(), TempCone) == Result.end() ) {
				Result.push_back(TempCone);
			};	
			ListAppendTime += double(clock() - begin);

			set<int>::iterator NeighborItr;
			begin = clock();
			for(NeighborItr=EdgeToTest.NeighborIndices.begin();NeighborItr!=EdgeToTest.NeighborIndices.end(); NeighborItr++){
				int Neighbor = *NeighborItr;
				if ( find(PretropGraphEdges.begin(), PretropGraphEdges.end(), Neighbor) == PretropGraphEdges.end() ) {
					if ( find(NotPretropGraphEdges.begin(), NotPretropGraphEdges.end(), Neighbor) == NotPretropGraphEdges.end() ) {
						if ( find(EdgesToTest.begin(), EdgesToTest.end(), Neighbor) == EdgesToTest.end() ) {
							EdgesToTest.push_back(Neighbor);
						};
					};
				};
			};
		} else {
			NotPretropGraphEdges.insert(EdgeToTestIndex);
		};
	};
	return Result;
}

//------------------------------------------------------------------------------				
int main(int argc, char* argv[]) {
	if (argc != 3) {
		cout << "Internal error: expected two arguments." << endl;
		return 1;
	}
	int n = atoi(argv[1]);
	bool Reduced = true;
	vector<Hull> Hulls;
	vector<vector<vector<GMP_Integer> > > Cyc4 = CyclicN(n,Reduced);
	if (Reduced == true) {
		n = n - 1;
	};
	vector<vector<vector<GMP_Integer> > >::iterator CycIt;
	int TestIndex = 0;
	for (CycIt=Cyc4.begin(); CycIt != Cyc4.end(); CycIt++) {
		Hulls.push_back(NewHull(*CycIt));
		TestIndex++;
	}
	vector<C_Polyhedron> ConeVector;	
	
	if (string(argv[2]) == "refinement") {
		// Start by initializing the objects.
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
					ConeIntersectionCount++;
					if ( find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
						NewCones.push_back(NewCone);
					};
				};
			};
			ConeVector = NewCones;
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", TreeLevel, Cones.size(), ConeVector.size(),ConeIntersectionCount);
			TreeLevel++;
		};
	}
	else if (string(argv[2]) == "new") {
		for (int HullIndex = 0; HullIndex != Hulls.size(); HullIndex++) {
			Hull H = Hulls[HullIndex];
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

					C_Polyhedron NewCone = *ConeIterator;
					
					vector<C_Polyhedron> ConesToAdd = TraverseEdgeSkeleton(NewCone, Pts, PointToIndexMap, Edges);

					for(int i = 0; i != ConesToAdd.size(); i++) {
						C_Polyhedron ConeToAdd = ConesToAdd[i];
						clock_t begin = clock();
						if ( find(NewCones.begin(), NewCones.end(), ConeToAdd) == NewCones.end() ) {
							NewCones.push_back(ConeToAdd);
						};
						ListAppendTime += clock() - begin;
					}
				};
				ConeVector = NewCones;
			};
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", HullIndex, Hulls.size(), ConeVector.size(),ConeIntersectionCount);
		};
		cout << "ConeVector count: " << ConeVector.size() << endl;
	};

	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;
	vector<Generator> gv;
	vector<C_Polyhedron>::iterator PolyItr;
	vector<vector<vector<GMP_Integer> > > InitialForms;

	// This is parsing and displaying all of the pretropisms.
	for (PolyItr=ConeVector.begin(); PolyItr != ConeVector.end(); PolyItr++) {
		Generator_System gs = (*PolyItr).minimized_generators();
		for (Generator_System::const_iterator gsi = gs.begin(),
		gs_end = gs.end(); gsi != gs_end; ++gsi) {
			Generator gen = *gsi;
			if ( find(gv.begin(), gv.end(), gen) == gv.end() ) {
				
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

					vector<vector<GMP_Integer> >::iterator InFormIter;
					for (InFormIter=InForm.begin(); InFormIter!=InForm.end(); InFormIter++) {
						PtIndices.push_back(MyHull.PointToIndexMap[*InFormIter]);
					}
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
	cout << "Intersection time: " << GetPolyhedralIntersectionTime() / CLOCKS_PER_SEC << endl;	
	cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
	cout << "List append time: " << ListAppendTime / CLOCKS_PER_SEC << endl;
}
