#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"
#include <algorithm>
#include "tbb/tbb.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <vector>

using namespace std;
using namespace Parma_Polyhedra_Library;
using namespace tbb;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

typedef vector<C_Polyhedron> ConeList;

double ListAppendTime;
double ContainmentTime;
double IntTime;
double OtherIntTime;
double streamtime;
int ConeIntersectionCount;

struct c_poly_compare {
    bool operator() (const C_Polyhedron& lhs, const C_Polyhedron& rhs) const{
        stringstream s1,s2;
        s1 << lhs.minimized_generators();
        s2 << rhs.minimized_generators();
        return s1.str() < s2.str();
    }
};

set<C_Polyhedron, c_poly_compare> ConeVector;
vector<vector<vector<GMP_Integer> > > CyclicN(int n, bool Reduced) {
	vector<vector<vector<GMP_Integer> > > System;
	int Length = n;
	if (Reduced == true) {
		Length = n - 1;
	};
	
	for (size_t i = 0; i != n-1; i++) {
		vector<vector<GMP_Integer> > Equation;
		for (size_t j = 0; j != n; j++) {
			vector<GMP_Integer> Monomial;
			set<GMP_Integer> OneSet;
			for (size_t k = 0; k != i+1; k++) {
				OneSet.insert((j+k)%n);
			};
			for (size_t l = 0; l != Length; l++) {
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


/*class IntersectTest {
	NewCone C_Polyhedron;
	Pts vector<vector<GMP_Integer> >;
	PointToIndexMap map<vector<GMP_Integer>,GMP_Integer>;
	Edges vector<Edge>;
	void operator() ( const blocked_range<size_t>& r ) const {
	
	
	}
}*/

//------------------------------------------------------------------------------
vector<set<C_Polyhedron,c_poly_compare> > IntersectCones(Hull H, set<C_Polyhedron,c_poly_compare> TestCones) {
	vector<vector<GMP_Integer> > Pts = H.Points;
	map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap = H.PointToIndexMap;
	vector<Edge> Edges = H.Edges;
	
	vector<set<C_Polyhedron,c_poly_compare> > Result;
	set<C_Polyhedron,c_poly_compare>::iterator it;
	for (it = TestCones.begin(); it != TestCones.end(); it++) {

		C_Polyhedron NewCone = *it;
		
		//take a random vector from cone
		vector<vector<GMP_Integer> > Rays = GeneratorSystemToPoints(NewCone.minimized_generators());
		vector<GMP_Integer> RandomVector(Rays[0].size(),0);
		for (size_t i = 0; i != Rays.size(); i++) {
			vector<GMP_Integer> Ray = Rays[i];
			for (size_t j = 0; j != Ray.size(); j++) {
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
		for(size_t EdgeIndex = 0; EdgeIndex != Edges.size(); EdgeIndex++) {
			if (SetsDoIntersect(InitialIndices, Edges[EdgeIndex].PointIndices)) {
				EdgesToTest.push_back(EdgeIndex);
			};
		};
		set<C_Polyhedron, c_poly_compare> NewCones;
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
				NewCones.insert(NewCone);
				break;
			};
			clock_t OtherIntBegin = clock();
			C_Polyhedron TempCone = IntersectCones(EdgeToTest.Cone, NewCone);
			OtherIntTime += double(clock() - OtherIntBegin);
			ConeIntersectionCount++;
			int ExpectedDimension = min(EdgeToTest.Cone.affine_dimension(),NewCone.affine_dimension()) - 1;
			if (ExpectedDimension <= TempCone.affine_dimension()) {
				PretropGraphEdges.insert(EdgeToTestIndex);
			
				clock_t begin = clock();
				NewCones.insert(TempCone);
				ListAppendTime += double(clock() - begin);

				set<int>::iterator NeighborItr;
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
		Result.push_back(NewCones);
	};
	return Result;
}

//------------------------------------------------------------------------------				
bool HullSort(Hull a, Hull b) {
    return a.Dimension > b.Dimension;
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
	

	//sort( Hulls.begin(), Hulls.end(), HullSort);
	
	if (string(argv[2]) == "refinement") {
		// Start by initializing the objects.
		vector<set<C_Polyhedron,c_poly_compare> > Cones;
		int HullIndex = 0;
		vector<Hull>::iterator it;
		for (it=Hulls.begin(); it != Hulls.end(); it++) {
			vector<Edge> Edges = (*it).Edges;
			set<C_Polyhedron,c_poly_compare> HullCones;
			vector<Edge>::iterator itr;

			for (itr=Edges.begin(); itr != Edges.end(); itr++) {
				HullCones.insert((*itr).Cone);
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
		vector<set<C_Polyhedron,c_poly_compare> >::iterator ConesItr;
		int TreeLevel = 1;
		for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
			set<C_Polyhedron,c_poly_compare> TestCones = *ConesItr;
			set<C_Polyhedron,c_poly_compare> NewCones;
			//Iterate through Cones
			set<C_Polyhedron,c_poly_compare>::iterator ConeItr;
			for (ConeItr=ConeVector.begin(); ConeItr != ConeVector.end(); ConeItr++) {
				//Iterate through TestCones
				C_Polyhedron Cone = *ConeItr;
				set<C_Polyhedron,c_poly_compare>::iterator TestConesItr;
				for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
					C_Polyhedron NewCone = IntersectCones(*TestConesItr, Cone);
					ConeIntersectionCount++;
					if ( find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
						NewCones.insert(NewCone);
					};
				};
			};
			ConeVector = NewCones;
			printf("Finished level %d of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", TreeLevel, Cones.size(), ConeVector.size(),ConeIntersectionCount);
			TreeLevel++;
		};
	}
	else if (string(argv[2]) == "new") {
		for (size_t HullIndex = 0; HullIndex != Hulls.size(); HullIndex++) {
			Hull H = Hulls[HullIndex];
			vector<Edge> Edges = H.Edges;
			vector<Facet> Facets = H.Facets;
			map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap = H.PointToIndexMap;
			vector<vector<GMP_Integer> > Pts = H.Points;
			
			//Initialize ConeVector if it's the first convex hull.
			if (HullIndex == 0) {
				vector<Edge>::iterator itr;
				for (itr=Edges.begin(); itr != Edges.end(); itr++) {
					ConeVector.insert((*itr).Cone);
				};
			} else {
				set<C_Polyhedron, c_poly_compare> NewCones;
				vector<C_Polyhedron>::iterator ConeIterator;
				


				vector<set<C_Polyhedron, c_poly_compare> >Arr;
				Arr = IntersectCones(H,ConeVector);
				for(size_t i = 0; i != ConeVector.size(); i++) {
					set<C_Polyhedron, c_poly_compare> ConesToAdd = Arr[i];
					set<C_Polyhedron, c_poly_compare>::iterator CPolyItr;
					for(CPolyItr = ConesToAdd.begin(); CPolyItr != ConesToAdd.end(); CPolyItr++) {
						C_Polyhedron ConeToAdd = *CPolyItr;
						clock_t begin = clock();
						NewCones.insert(ConeToAdd);
						ListAppendTime += clock() - begin;
					}
				};
				vector<Constraint> Constraints;/*
				for(size_t i = 0; i != NewCones.size(); i++){
					Constraint_System cstest = NewCones[i].minimized_constraints();
					for (Constraint_System::const_iterator it = cstest.begin(),
					cs_end = cstest.end(); it != cs_end; ++it) {
						if ( find(Constraints.begin(), Constraints.end(), *it) == Constraints.end() ) {
							Constraints.push_back(*it);
						};
					};
				};*/
				cout << "Intersection time: " << GetPolyhedralIntersectionTime() / CLOCKS_PER_SEC << endl;	
				cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
				cout << "List append time: " << ListAppendTime / CLOCKS_PER_SEC << endl;
				cout << "ConstraintCount: " << Constraints.size() << endl;
				ConeVector = NewCones;
			};
			printf("Finished level %lu of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", HullIndex, Hulls.size(), ConeVector.size(),ConeIntersectionCount);
		};
		cout << "ConeVector count: " << ConeVector.size() << endl;
	};

	cout << "Finished intersecting C_Polyhedrons." << endl << endl << endl;
	vector<Generator> gv;
	set<C_Polyhedron,c_poly_compare>::iterator PolyItr;
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
	cout << "IntTime: " << IntTime / CLOCKS_PER_SEC << endl;
	cout << "IntTime: " << OtherIntTime / CLOCKS_PER_SEC << endl;
	cout << "Stream time: " << streamtime / CLOCKS_PER_SEC << endl;
}
