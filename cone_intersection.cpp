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

double ListAppendTime, ContainmentTime, IntersectionTime, ParallelTime;
int ConeIntersectionCount;
vector<Cone> ConeVector;

//------------------------------------------------------------------------------
struct c_poly_compare {
    bool operator() (const C_Polyhedron& lhs, const C_Polyhedron& rhs) const{
        stringstream s1,s2;
        s1 << lhs.minimized_generators();
        s2 << rhs.minimized_generators();
        return s1.str() < s2.str();
    }
};

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------				
bool ConeSort(C_Polyhedron a, C_Polyhedron b) {
	return a.affine_dimension() > b.affine_dimension();
}

// Globally:
bool operator == (const Cone &L, const Cone &R)
{
    return(L.Polyhedron == R.Polyhedron);
}

//------------------------------------------------------------------------------
class IntersectTest {
	int SpaceDimension;
	int HullIndex;
	vector<Edge> &Edges;
	vector<vector<GMP_Integer> > &Pts;
	map<vector<GMP_Integer>,GMP_Integer> &PointToIndexMap;
	vector<Cone> &TestCones;
	vector<vector<Cone> > &Result;
	public:
		IntersectTest(int dim, int hi, vector<Edge> &e, vector<vector<GMP_Integer> > &pts, map<vector<GMP_Integer>,GMP_Integer> &ptmap, vector<Cone> &cones, vector<vector<Cone> > &output): SpaceDimension(dim), HullIndex(hi), Edges(e), Pts(pts), PointToIndexMap(ptmap), TestCones(cones), Result(output) { }
    void operator() (const blocked_range<size_t>& r ) const	
{
	for(size_t ConeIndex=r.begin(); ConeIndex != r.end(); ConeIndex++) {
	
		Cone NewCone = TestCones[ConeIndex];
		int NumberOfHulls = NewCone.IntersectionIndices.size();
		set<int> NewConeIntersectionIndices = NewCone.IntersectionIndices[HullIndex];
		//take a random vector from cone		
		vector<GMP_Integer> RandomVector(SpaceDimension, 0);
		Generator_System gs = NewCone.Polyhedron.minimized_generators();
		for (Generator_System::const_iterator i = gs.begin(),
		gs_end = gs.end(); i != gs_end; ++i) {
			for (size_t j = 0; j != SpaceDimension; j++) {
				RandomVector[j] += (*i).coefficient(Variable(j));
			};
		};
		//take initial form using that vector.
		vector<vector<GMP_Integer> > InitialForm = FindInitialForm(Pts, RandomVector);
	
		set<GMP_Integer> InitialIndices;
		for (vector<vector<GMP_Integer> >::iterator InitialFormItr=InitialForm.begin();
		InitialFormItr != InitialForm.end(); InitialFormItr++) {
			InitialIndices.insert(PointToIndexMap[*InitialFormItr]);
		};
	
		// List of indices of edges to visit.
		vector<int> InitialEdgesToTest;
		for(size_t EdgeIndex = 0; EdgeIndex != Edges.size(); EdgeIndex++) {
			if (SetsDoIntersect(InitialIndices, Edges[EdgeIndex].PointIndices)) {
				InitialEdgesToTest.push_back(EdgeIndex);
			};
		};
		
		vector<Cone> NewCones;
		// Explore edge skeleton
		set<int> PretropGraphEdges;
		set<int> NotPretropGraphEdges;
		vector<int> EdgesToTest;
		while(!EdgesToTest.empty() || !InitialEdgesToTest.empty()) {
			int EdgeToTestIndex;
			bool DoConeContainment;
			if (!InitialEdgesToTest.empty()) {
				EdgeToTestIndex = InitialEdgesToTest.back();
				InitialEdgesToTest.pop_back();
				DoConeContainment = true;
			} else {
				EdgeToTestIndex = EdgesToTest.back();
				EdgesToTest.pop_back();
				DoConeContainment = false;
			};
			Edge EdgeToTest = Edges[EdgeToTestIndex];
			
			clock_t ContainmentBegin = clock();
			// Only do the containment check if it's one of the initial edges.
			bool ConeDoesContain = false;
			if (DoConeContainment) {
				ConeDoesContain = EdgeToTest.EdgeCone.Polyhedron.contains(NewCone.Polyhedron);
			};
			ContainmentTime += double(clock() - ContainmentBegin);

			if (ConeDoesContain) {
				if (find(NewCones.begin(), NewCones.end(), NewCone) == NewCones.end() ) {
					NewCones.push_back(NewCone);
				};
				break;
			};
			clock_t IntBegin = clock();

			Constraint_System cs1 = EdgeToTest.EdgeCone.Polyhedron.minimized_constraints();
			Constraint_System cs2 = NewCone.Polyhedron.minimized_constraints();
			for (Constraint_System::const_iterator i = cs1.begin(),
			cs1_end = cs1.end(); i != cs1_end; ++i) {
				cs2.insert(*i);
			};
			Recycle_Input dummy;
			Cone TempCone;
			C_Polyhedron TempCPolyhedron(cs2, dummy);
			TempCone.Polyhedron = TempCPolyhedron;
			TempCone.Polyhedron.affine_dimension();
			IntersectionTime += double(clock() - IntBegin);
			ConeIntersectionCount++;
			int ExpectedDimension = min(EdgeToTest.EdgeCone.Polyhedron.affine_dimension(),NewCone.Polyhedron.affine_dimension()) - 1;
			if (ExpectedDimension <= TempCone.Polyhedron.affine_dimension()) {
				PretropGraphEdges.insert(EdgeToTestIndex);
			
				clock_t begin = clock();
				vector<Cone>:: iterator ConeIt;
				
				ConeIt = find(NewCones.begin(), NewCones.end(), TempCone);
				if (ConeIt == NewCones.end()) {
					// Initialize the sets by intersecting the sets of EdgeToTest.EdgeCone and NewCone
					vector<set<int> > InitialSet (NewCone.IntersectionIndices.size());
					TempCone.IntersectionIndices = InitialSet;
					for (size_t i = HullIndex + 1; i != NewCone.IntersectionIndices.size(); i++) {
						TempCone.IntersectionIndices[i] = IntersectSets(EdgeToTest.EdgeCone.IntersectionIndices[i], NewCone.IntersectionIndices[i]);
					};
					NewCones.push_back(TempCone);
				} else {
					for (size_t i = HullIndex + 1; i != NewCone.IntersectionIndices.size(); i++) {
						(*ConeIt).IntersectionIndices[i] = IntersectSets((*ConeIt).IntersectionIndices[i], NewCone.IntersectionIndices[i]);
						(*ConeIt).IntersectionIndices[i] = IntersectSets((*ConeIt).IntersectionIndices[i], EdgeToTest.EdgeCone.IntersectionIndices[i]);
					};

				};
				
				
				ListAppendTime += double(clock() - begin);

				set<int>::iterator NeighborItr;
				for(NeighborItr=EdgeToTest.NeighborIndices.begin();NeighborItr!=EdgeToTest.NeighborIndices.end(); NeighborItr++){
					int Neighbor = *NeighborItr;
					if (( find(NewConeIntersectionIndices.begin(), NewConeIntersectionIndices.end(), Neighbor) != NewConeIntersectionIndices.end())
					&& ( find(PretropGraphEdges.begin(), PretropGraphEdges.end(), Neighbor) == PretropGraphEdges.end() )
					&& ( find(NotPretropGraphEdges.begin(), NotPretropGraphEdges.end(), Neighbor) == NotPretropGraphEdges.end() )
					&& ( find(EdgesToTest.begin(), EdgesToTest.end(), Neighbor) == EdgesToTest.end() )
					&& ( find(InitialEdgesToTest.begin(), InitialEdgesToTest.end(), Neighbor) == InitialEdgesToTest.end() )) {
						EdgesToTest.push_back(Neighbor);
					};
				};
			} else {
				NotPretropGraphEdges.insert(EdgeToTestIndex);
			};
		};
		Result[ConeIndex] = NewCones;
	};
}	
};

//------------------------------------------------------------------------------
bool HullSort(Hull a, Hull b) {
	return a.AffineDimension > b.AffineDimension;
}

//------------------------------------------------------------------------------				
int main(int argc, char* argv[]) {
	clock_t StartTime = clock();
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
		vector<vector<Cone> > Cones;
		int HullIndex = 0;
		vector<Hull>::iterator it;
		for (it=Hulls.begin(); it != Hulls.end(); it++) {
			vector<Edge> Edges = (*it).Edges;
			vector<Cone> HullCones;
			vector<Edge>::iterator itr;

			for (itr=Edges.begin(); itr != Edges.end(); itr++) {
				HullCones.push_back((*itr).EdgeCone);
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
		vector<vector<Cone> >::iterator ConesItr;
		int TreeLevel = 1;
		for (ConesItr=Cones.begin(); ConesItr != Cones.end(); ConesItr++) {
			vector<Cone> TestCones = *ConesItr;
			vector<Cone> NewCones;
			//Iterate through Cones
			vector<Cone>::iterator ConeItr;
			for (ConeItr=ConeVector.begin(); ConeItr != ConeVector.end(); ConeItr++) {
				//Iterate through TestCones
				Cone Cone1 = *ConeItr;
				vector<Cone>::iterator TestConesItr;
				for (TestConesItr=TestCones.begin(); TestConesItr != TestCones.end(); TestConesItr++) {
					Cone NewCone = IntersectCones(*TestConesItr, Cone1);
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
		for(int i = 0; i != Hulls.size(); i++){
			for(int j = 0; j != Hulls[i].Edges.size(); j++){
				vector<set<int> > InitialSet (Hulls.size());
				Hulls[i].Edges[j].EdgeCone.IntersectionIndices = InitialSet;
			};
		};	
		clock_t PreintTimeStart = clock();
		int TotalInt = 0;
		int NonInt = 0;
		int ExpectedDim = Hulls[0].Edges[0].EdgeCone.Polyhedron.affine_dimension() - 1;
		// TODO: switch this to explore the edge skeleton method.
		for(int i = 0; i != Hulls.size(); i++){
			vector<Edge> Edges1 = Hulls[i].Edges;
			for(int j = i+1; j != Hulls.size(); j++){
				vector<Edge> Edges2 = Hulls[j].Edges;
			
				for(int k = 0; k != Edges1.size(); k++){
					for(int l = 0; l != Edges2.size(); l++){
						bool ConesDoIntersect = IntersectCones(Edges1[k].EdgeCone.Polyhedron, Edges2[l].EdgeCone.Polyhedron).affine_dimension() >= ExpectedDim;
						if (ConesDoIntersect) {
							Hulls[i].Edges[k].EdgeCone.IntersectionIndices[j].insert(l);
							Hulls[j].Edges[l].EdgeCone.IntersectionIndices[i].insert(k);
						};
						if (!ConesDoIntersect) {
							NonInt++;
						};
						TotalInt++;
					};
				};
			};
			printf("Finished level %d of pre-intersections.\n", i);
		};
		cout << "Total Intersections: " << TotalInt << ", Non Intersections: " << NonInt << endl;
		cout << "Preintersection time: " << double(clock() - PreintTimeStart) / CLOCKS_PER_SEC << endl;

	
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
					ConeVector.push_back((*itr).EdgeCone);
				};
			} else {
				vector<Cone> NewCones;
				vector<Cone>::iterator ConeIterator;
	
	
				vector<vector<Cone> > Arr(ConeVector.size());
				task_scheduler_init init(1);

				clock_t ParallelBegin = clock();
				parallel_for(blocked_range<size_t>(0,ConeVector.size()),IntersectTest(H.SpaceDimension, HullIndex,Edges,Pts,PointToIndexMap,ConeVector,Arr));
				ParallelTime += double(clock() - ParallelBegin);
				for(size_t i = 0; i != ConeVector.size(); i++) {
					vector<Cone> ConesToAdd = Arr[i];
					vector<Cone>::iterator CPolyItr;
					for(CPolyItr = ConesToAdd.begin(); CPolyItr != ConesToAdd.end(); CPolyItr++) {
						Cone ConeToAdd = *CPolyItr;
						clock_t begin = clock();

						vector<Cone>:: iterator ConeIt = find(NewCones.begin(), NewCones.end(), ConeToAdd);
						if (ConeIt == NewCones.end()) {
							NewCones.push_back(ConeToAdd);
						} else {
							for (size_t i = HullIndex + 1; i != ConeToAdd.IntersectionIndices.size(); i++) {
								(*ConeIt).IntersectionIndices[i] = IntersectSets((*ConeIt).IntersectionIndices[i], ConeToAdd.IntersectionIndices[i]);
							};
						};
						ListAppendTime += clock() - begin;
					};
				};
				vector<Constraint> Constraints;
				/*
				//TODO: constraint relation table.
				//TODO: think about using affine_image to do dimension reduction.
				//TODO: try using standard integer instead of GMP_Integer
				for(size_t i = 0; i != NewCones.size(); i++){
					Constraint_System cstest = NewCones[i].minimized_constraints();
					for (Constraint_System::const_iterator it = cstest.begin(),
					cs_end = cstest.end(); it != cs_end; ++it) {
						if ( find(Constraints.begin(), Constraints.end(), *it) == Constraints.end() ) {
							Constraints.push_back(*it);
						};
					};
				};*/
				cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
				cout << "List append time: " << ListAppendTime / CLOCKS_PER_SEC << endl;
				cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
				cout << "Parallel time: " << ParallelTime / CLOCKS_PER_SEC << endl;
				cout << "Total elapsed time: " << double(clock() - StartTime) / CLOCKS_PER_SEC << endl;
				//cout << "ConstraintCount: " << Constraints.size() << endl;
				ConeVector = NewCones;
			};
			printf("Finished level %lu of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", HullIndex, Hulls.size(), ConeVector.size(),ConeIntersectionCount);
		};
		cout << "ConeVector count: " << ConeVector.size() << endl;
	}
	else if (string(argv[2]) == "graph") {
		for(int i = 0; i != Hulls.size(); i++){
			for(int j = 0; j != Hulls[i].Edges.size(); j++){
				vector<set<int> > InitialSet (Hulls.size());
				Hulls[i].Edges[j].EdgeCone.IntersectionIndices = InitialSet;
			};
		};	
		clock_t PreintTimeStart = clock();
		int TotalInt = 0;
		int NonInt = 0;
		int ExpectedDim = Hulls[0].Edges[0].EdgeCone.Polyhedron.affine_dimension() - 1;
		// TODO: switch this to explore the edge skeleton method.
		for(int i = 0; i != Hulls.size(); i++){
			vector<Edge> Edges1 = Hulls[i].Edges;
			for(int j = i+1; j != Hulls.size(); j++){
				vector<Edge> Edges2 = Hulls[j].Edges;
			
				for(int k = 0; k != Edges1.size(); k++){
					for(int l = 0; l != Edges2.size(); l++){
						bool ConesDoIntersect = IntersectCones(Edges1[k].EdgeCone.Polyhedron, Edges2[l].EdgeCone.Polyhedron).affine_dimension() >= ExpectedDim;
						if (ConesDoIntersect) {
							Hulls[i].Edges[k].EdgeCone.IntersectionIndices[j].insert(l);
							Hulls[j].Edges[l].EdgeCone.IntersectionIndices[i].insert(k);
						};
						if (!ConesDoIntersect) {
							NonInt++;
						};
						TotalInt++;
					};
				};
			};
			printf("Finished level %d of pre-intersections.\n", i);
		};
		cout << "Total Intersections: " << TotalInt << ", Non Intersections: " << NonInt << endl;
		cout << "Preintersection time: " << double(clock() - PreintTimeStart) / CLOCKS_PER_SEC << endl;
		vector<vector<int> > Cliques;
		for (size_t i = 0; i != Hulls[0].Edges.size(); i++) {
			vector<int> Clique;
			Clique.push_back(i);
			Cliques.push_back(Clique);
		};
		
		for (size_t i = 1; i != Hulls.size(); i++) {
			vector<vector<int> > NewCliques;
			for(size_t j = 0; j != Cliques.size(); j++) {
				vector<int> Clique = Cliques[j];
				for(size_t k = 0; k != Hulls[i].Edges.size(); k++) {
					Edge E = Hulls[i].Edges[k];
					bool IsPartOfClique = true;
					for (size_t m = 0; m != Clique.size(); m++) {
						set<int> CliqueSet = E.EdgeCone.IntersectionIndices[m];
						if (find(CliqueSet.begin(), CliqueSet.end(), Clique[m]) == CliqueSet.end()) {
							IsPartOfClique = false;
							break;
						};
					};
					if (IsPartOfClique) {
						vector<int> NewClique (Clique.size());
						for (int l = 0; l != Clique.size(); l++) {
							NewClique[l] = Clique[l];
						};
						NewClique.push_back(k);
						NewCliques.push_back(NewClique);
					};
				};
			};
			cout << "NewCliques count: " << NewCliques.size() << endl;
			Cliques = NewCliques;
		};
		int NonIntersection = 0;
		for (size_t i = 0; i != Cliques.size(); i++) {
			vector<int> Clique = Cliques[i];
			C_Polyhedron ph = Hulls[0].Edges[Clique[0]].EdgeCone.Polyhedron;
			bool ShouldAddCone = true;
			for (size_t j = 1; j != Clique.size(); j++) {
				if (ph.affine_dimension() == 0) {
					ShouldAddCone = false;
					break;
				}
				ph = IntersectCones(ph,Hulls[j].Edges[Clique[j]].EdgeCone.Polyhedron);
			};
			if (ph.affine_dimension() == 0) {
				NonIntersection++;
			};
			if (ShouldAddCone) {
				Cone NewCone;
				NewCone.Polyhedron = ph;
				ConeVector.push_back(NewCone);
			};
		};
		cout << "NonIntersection: " << NonIntersection << endl;
	};
	
	cout << "Finished intersecting cones." << endl << endl << endl;
	vector<Generator> gv;
	vector<Cone>::iterator PolyItr;
	vector<vector<vector<GMP_Integer> > > InitialForms;

	// This is parsing and displaying all of the pretropisms.
	for (PolyItr=ConeVector.begin(); PolyItr != ConeVector.end(); PolyItr++) {
		Generator_System gs = (*PolyItr).Polyhedron.minimized_generators();
		for (Generator_System::const_iterator gsi = gs.begin(),
		gs_end = gs.end(); gsi != gs_end; ++gsi) {
			Generator gen = *gsi;
			if (gen.is_point() or gen.is_line()) {
				continue;
			};
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
//					cout << "Pretropism: " << gen << endl;
				}
			};
		};
	};
	cout << "Number of pretropisms found: " << gv.size() << endl;
	cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
	cout << "List append time: " << ListAppendTime / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Parallel time: " << ParallelTime / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;
}
