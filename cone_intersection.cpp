#include <iostream>
#include <ppl.hh>
#include <ctime>
#include "prevariety_util.h"
#include "cone_tree.h"
#include <algorithm>
#include "tbb/tbb.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <vector>
#include <string>
#include <sstream>

using namespace std;
using namespace Parma_Polyhedra_Library;
using namespace tbb;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

double TreeTime1, TreeTime2, ContainmentTime, IntersectionTime, ParallelTime, TestTime, GetConesTime;
int ConeIntersectionCount;
vector<Cone> ConeVector;
//------------------------------------------------------------------------------				
void PreintersectCones (Hull &H1, Hull &H2) {

};

//------------------------------------------------------------------------------
struct c_poly_compare {
    bool operator() (const NNC_Polyhedron& lhs, const NNC_Polyhedron& rhs) const{
        stringstream s1,s2;
        s1 << lhs.minimized_generators();
        s2 << rhs.minimized_generators();
        return s1.str() < s2.str();
    }
};

//------------------------------------------------------------------------------
vector<vector<vector<int> > > CyclicN(int n, bool Reduced) {
	vector<vector<vector<int> > > System;
	int Length = n;
	if (Reduced == true) {
		Length = n - 1;
	};
	
	for (size_t i = 0; i != n-1; i++) {
		vector<vector<int> > Equation;
		for (size_t j = 0; j != n; j++) {
			vector<int> Monomial;
			set<int> OneSet;
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
bool ConeSort(NNC_Polyhedron a, NNC_Polyhedron b) {
	return a.affine_dimension() > b.affine_dimension();
}

//------------------------------------------------------------------------------
bool operator == (const Cone &L, const Cone &R)
{
    return(L.Polyhedron == R.Polyhedron);
}
/*
//------------------------------------------------------------------------------
class IntersectTest {
	int SpaceDimension;
	int HullIndex;
	vector<Edge> &Edges;
	vector<vector<int> > &Pts;
	map<vector<int>,int> &PointToIndexMap;
	vector<Cone> &TestCones;
	vector<vector<Cone> > &Result;
	public:
		IntersectTest(int dim, int hi, vector<Edge> &e, vector<vector<int> > &pts, map<vector<int>,int> &ptmap, vector<Cone> &cones, vector<vector<Cone> > &output): SpaceDimension(dim), HullIndex(hi), Edges(e), Pts(pts), PointToIndexMap(ptmap), TestCones(cones), Result(output) { }
    void operator() (const blocked_range<size_t>& r ) const	
{
	for(size_t ConeIndex=r.begin(); ConeIndex != r.end(); ConeIndex++) {
		Recycle_Input dummy;
		Node Tree;
		Tree.IsLeaf = false;
		Cone NewCone = TestCones[ConeIndex];
		set<int> NewConeIntersectionIndices = NewCone.IntersectionIndices[HullIndex];
		//take a random vector from cone		
		vector<int> RandomVector(SpaceDimension, 0);
		Generator_System gs = NewCone.Polyhedron.minimized_generators();
		for (Generator_System::const_iterator i = gs.begin(),
		gs_end = gs.end(); i != gs_end; ++i) {
			for (size_t j = 0; j != SpaceDimension; j++) {
				stringstream s;
				s << (*i).coefficient(Variable(j));
				int ToAppend;
				istringstream(s.str()) >> ToAppend;
				RandomVector[j] += ToAppend;
			};
		};
		//take initial form using that vector.
		vector<vector<int> > InitialForm = FindInitialForm(Pts, RandomVector);
	
		set<int> InitialIndices;
		for (vector<vector<int> >::iterator InitialFormItr=InitialForm.begin();
		InitialFormItr != InitialForm.end(); InitialFormItr++) {
			InitialIndices.insert(PointToIndexMap[*InitialFormItr]);
		};
	
		// List of indices of edges to visit.
		vector<int> InitialEdgesToTest;
		for(size_t EdgeIndex = 0; EdgeIndex != Edges.size(); EdgeIndex++) {
			if (SetsDoIntersect(InitialIndices, Edges[EdgeIndex].PointIndices)
			&& ( find(NewConeIntersectionIndices.begin(), NewConeIntersectionIndices.end(), EdgeIndex) != NewConeIntersectionIndices.end())) {
				InitialEdgesToTest.push_back(EdgeIndex);
			};
		};
		
		vector<Cone> NewCones;
		// Explore edge skeleton
		set<int> PretropGraphEdges;
		set<int> NotPretropGraphEdges;
		vector<int> EdgesToTest;
		Edge EdgeToTest;
		int EdgeToTestIndex;
		bool DoConeContainment;
		while(!EdgesToTest.empty() || !InitialEdgesToTest.empty()) {
			if (!InitialEdgesToTest.empty()) {
				EdgeToTestIndex = InitialEdgesToTest.back();
				InitialEdgesToTest.pop_back();
				DoConeContainment = true;
			} else {
				EdgeToTestIndex = EdgesToTest.back();
				EdgesToTest.pop_back();
				DoConeContainment = false;
			};
			clock_t BeginTime = clock();
			EdgeToTest = Edges[EdgeToTestIndex];
			TestTime += double(clock() - BeginTime);
			clock_t ContainmentBegin = clock();
			// Only do the containment check if it's one of the initial edges.
			bool ConeDoesContain = false;
			if (DoConeContainment) {
				ConeDoesContain = EdgeToTest.EdgeCone.Polyhedron.contains(NewCone.Polyhedron);
			};
			ContainmentTime += double(clock() - ContainmentBegin);
			
			if (ConeDoesContain) {
				clock_t begin = clock();
				InsertCone(Tree,NewCone,HullIndex+1);
				TreeTime1 += double(clock() - begin);
				break;
			};
			
			clock_t IntBegin = clock();

			Constraint_System cs1 = EdgeToTest.EdgeCone.Polyhedron.minimized_constraints();
			Constraint_System cs2 = NewCone.Polyhedron.minimized_constraints();
			
			for (Constraint_System::const_iterator i = cs1.begin(),
			cs1_end = cs1.end(); i != cs1_end; ++i) {
				cs2.insert(*i);
			};
			Cone TempCone;
			NNC_Polyhedron TempCPolyhedron(cs2, dummy);
			TempCone.Polyhedron = TempCPolyhedron;
			TempCone.Polyhedron.affine_dimension();
			ConeIntersectionCount++;
			int ExpectedDimension = min(EdgeToTest.EdgeCone.Polyhedron.affine_dimension(),NewCone.Polyhedron.affine_dimension()) - 1;
			IntersectionTime += double(clock() - IntBegin);

			if (ExpectedDimension <= TempCone.Polyhedron.affine_dimension()) {
				PretropGraphEdges.insert(EdgeToTestIndex);
			
				clock_t begin = clock();
				vector<set<int> > InitialSet (NewCone.IntersectionIndices.size());
				TempCone.IntersectionIndices = InitialSet;
				for (size_t i = HullIndex + 1; i != NewCone.IntersectionIndices.size(); i++) {
					TempCone.IntersectionIndices[i] =  IntersectSets(EdgeToTest.EdgeCone.IntersectionIndices[i], NewCone.IntersectionIndices[i]);
				};
				InsertCone(Tree, TempCone,HullIndex+1);
				TreeTime1 += double(clock() - begin);
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

		Result[ConeIndex] = GetCones(Tree, HullIndex+1);
	};
}	
};*/

vector<Cone> IntersectTestTwo(int SpaceDimension, int HullIndex, vector<Edge> &Edges, vector<vector<int> > &Pts, map<vector<int>,int> &PointToIndexMap, vector<Cone> &TestCones, map<int,vector<int> > &IndexToPointMap)
{
	Node Tree;
	Tree.IsLeaf = false;
	for(size_t ConeIndex= 0; ConeIndex != TestCones.size(); ConeIndex++) {
		Recycle_Input dummy;
		Cone NewCone = TestCones[ConeIndex];
		set<int> NewConeIntersectionIndices = NewCone.IntersectionIndices[HullIndex];
		//take a random vector from cone		
		vector<int> RandomVector(SpaceDimension, 0);
		Generator_System gs = NewCone.Polyhedron.minimized_generators();
		for (Generator_System::const_iterator i = gs.begin(),gs_end = gs.end(); i != gs_end; ++i) {
			if (!(*i).is_ray() && !(*i).is_line()) {
				continue;
			};

			for (size_t j = 0; j != SpaceDimension; j++) {
				stringstream s;
				s << (*i).coefficient(Variable(j));
				int ToAppend;
				istringstream(s.str()) >> ToAppend;
				RandomVector[j] += ToAppend;
			};
		};
		//take initial form using that vector.
		vector<vector<int> > InitialForm = FindInitialForm(Pts, RandomVector);

		set<int> InitialIndices;
		for (vector<vector<int> >::iterator InitialFormItr=InitialForm.begin();
		InitialFormItr != InitialForm.end(); InitialFormItr++) {
			InitialIndices.insert(PointToIndexMap[*InitialFormItr]);
		};
	
		// List of indices of edges to visit.
		vector<int> InitialEdgesToTest;
		for(size_t EdgeIndex = 0; EdgeIndex != Edges.size(); EdgeIndex++) {
			if (SetsDoIntersect(InitialIndices, Edges[EdgeIndex].PointIndices)
			&& ( find(NewConeIntersectionIndices.begin(), NewConeIntersectionIndices.end(), EdgeIndex) != NewConeIntersectionIndices.end())
			) {
				InitialEdgesToTest.push_back(EdgeIndex);
			};
		};
		
		vector<Cone> NewCones;
		// Explore edge skeleton
		set<int> PretropGraphEdges;
		set<int> NotPretropGraphEdges;
		vector<int> EdgesToTest;
		Edge EdgeToTest;
		int EdgeToTestIndex;
		bool DoConeContainment;
		while(!EdgesToTest.empty() || !InitialEdgesToTest.empty()) {
			if (!InitialEdgesToTest.empty()) {
				EdgeToTestIndex = InitialEdgesToTest.back();
				InitialEdgesToTest.pop_back();
				DoConeContainment = true;
			} else {
				EdgeToTestIndex = EdgesToTest.back();
				EdgesToTest.pop_back();
				DoConeContainment = false;
			};
			clock_t BeginTime = clock();
			EdgeToTest = Edges[EdgeToTestIndex];
			TestTime += double(clock() - BeginTime);
			clock_t ContainmentBegin = clock();
			// Only do the containment check if it's one of the initial edges.
			bool ConeDoesContain = false;
			if (DoConeContainment) {
				ConeDoesContain = EdgeToTest.EdgeCone.Polyhedron.contains(NewCone.Polyhedron);
				ConeDoesContain = false;
			};
			ContainmentTime += double(clock() - ContainmentBegin);
			
			if (ConeDoesContain) {
				clock_t begin = clock();
				InsertCone(Tree,NewCone,HullIndex+1);
				TreeTime1 += double(clock() - begin);
				break;
			};
			
			clock_t IntBegin = clock();

			Constraint_System cs1 = EdgeToTest.EdgeCone.Polyhedron.minimized_constraints();
			Constraint_System cs2 = NewCone.Polyhedron.minimized_constraints();
			
			for (Constraint_System::const_iterator i = cs1.begin(),
			cs1_end = cs1.end(); i != cs1_end; ++i) {
				cs2.insert(*i);
			};
			Cone TempCone;
			NNC_Polyhedron TempCPolyhedron(cs2, dummy);
			TempCone.Polyhedron = TempCPolyhedron;
			ConeIntersectionCount++;
			IntersectionTime += double(clock() - IntBegin);
			if (TempCone.Polyhedron.affine_dimension() > 0) {
				PretropGraphEdges.insert(EdgeToTestIndex);
			
				clock_t begin = clock();
				vector<set<int> > InitialSet (NewCone.IntersectionIndices.size());
				TempCone.IntersectionIndices = InitialSet;
				for (size_t i = HullIndex + 1; i != NewCone.IntersectionIndices.size(); i++) {
					TempCone.IntersectionIndices[i] = IntersectSets(EdgeToTest.EdgeCone.IntersectionIndices[i], NewCone.IntersectionIndices[i]);
				};
				InsertCone(Tree, TempCone,HullIndex+1);
				TreeTime1 += double(clock() - begin);
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
	};
	clock_t BeginGetConesTime = clock();
	vector<Cone> Cones = GetCones(Tree, HullIndex + 1);
	GetConesTime += double(clock() - BeginGetConesTime);
	return Cones;
}


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
	vector<vector<vector<int> > > Cyc4 = CyclicN(n,Reduced);
	if (Reduced == true) {
		n = n - 1;
	};
	
	for (size_t i = 0; i != Cyc4.size(); i++) {
		if ((i != 0) and (i != 2)) {
			//continue;
		};
		if ((i != 0) and (i != 1) and (i != 2) and (i != 4)) {
			//continue;
		};
		bool UseHalfOpenCones;
		if (i==0) {
			UseHalfOpenCones = true;
		} else {
			UseHalfOpenCones = false;
		}
		Hulls.push_back(NewHull(Cyc4[i], UseHalfOpenCones));
	}
	//sort( Hulls.begin(), Hulls.end(), HullSort);
	
	double HullTime = double(clock() - StartTime);
	clock_t AlgorithmStartTime = clock();
	double PreintersectTime = 0;
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
	}	else if (string(argv[2]) == "new") {
		for(int i = 0; i != Hulls.size(); i++){
			for(int j = 0; j != Hulls[i].Edges.size(); j++){
				vector<set<int> > InitialSet (Hulls.size());
				Hulls[i].Edges[j].EdgeCone.IntersectionIndices = InitialSet;
			};
		};
		clock_t PreintTimeStart = clock();
		int TotalInt = 0;
		int NonInt = 0;
		int ExpectedDim = 1;
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
		PreintersectTime = double(clock() - PreintTimeStart);
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;

	
		for (size_t HullIndex = 0; HullIndex != Hulls.size(); HullIndex++) {
			Hull H = Hulls[HullIndex];
			vector<Edge> Edges = H.Edges;
			vector<Facet> Facets = H.Facets;
			map<vector<int>,int> PointToIndexMap = H.PointToIndexMap;
			vector<vector<int> > Pts = H.Points;
			
			//Initialize ConeVector if it's the first convex hull.
			if (HullIndex == 0) {
				vector<Edge>::iterator itr;
				for (itr=Edges.begin(); itr != Edges.end(); itr++) {
					ConeVector.push_back((*itr).EdgeCone);
				};
			} else {

				clock_t ParallelBegin = clock();
				vector<Cone> NewCones = IntersectTestTwo(H.SpaceDimension, HullIndex, Edges, Pts, PointToIndexMap, ConeVector, H.IndexToPointMap);
				
				ParallelTime += double(clock() - ParallelBegin);
				clock_t begin = clock();
				TreeTime2 += clock() - begin;
				clock_t GetConeTimeBegin = clock();
				GetConesTime += double(clock() - GetConeTimeBegin);
				vector<Constraint> Constraints;/*
				for(size_t i = 0; i != NewCones.size(); i++) {
					Constraint_System cstest = NewCones[i].Polyhedron.minimized_constraints();
					int count = 0;
					for (Constraint_System::const_iterator it = cstest.begin(),
					cs_end = cstest.end(); it != cs_end; ++it) {
						count++;
					};
					cout << count << endl;
				};
				cout << "Pause" << endl;
				cin.get();*/
				/*
				//TODO: constraint relation table.
				//TODO: think about using affine_image to do dimension reduction.
				//TODO: try using standard integer instead of int
				
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
				cout << "List append time for loop: " << TreeTime1 / CLOCKS_PER_SEC << endl;
				cout << "List append time not in for loop: " << TreeTime2 / CLOCKS_PER_SEC << endl;
				cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
				cout << "Parallel time: " << ParallelTime / CLOCKS_PER_SEC << endl;
				cout << "Total elapsed time: " << double(clock() - StartTime) / CLOCKS_PER_SEC << endl;
				//cout << "ConstraintCount: " << Constraints.size() << endl;
				ConeVector = NewCones;
			};
			printf("Finished level %lu of tree with %lu levels. %lu cones remain at this level. IntersectionCount = %d.\n", HullIndex, Hulls.size()-1, ConeVector.size(),ConeIntersectionCount);
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
		PreintersectTime = double(clock() - PreintTimeStart);
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
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
			NNC_Polyhedron ph = Hulls[0].Edges[Clique[0]].EdgeCone.Polyhedron;
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
	clock_t CleanupStart = clock();
	vector<Generator> gv;
	vector<Cone>::iterator PolyItr;
	double AlgTime = double(clock() - AlgorithmStartTime);
	// This is parsing and displaying all of the pretropisms.
	vector<vector<int> > Pretropisms;
	for (PolyItr=ConeVector.begin(); PolyItr != ConeVector.end(); PolyItr++) {
		Generator_System gs = (*PolyItr).Polyhedron.minimized_generators();
		for (Generator_System::const_iterator gsi = gs.begin(),
		gs_end = gs.end(); gsi != gs_end; ++gsi) {
			Generator gen = *gsi;
			if (gen.is_point() or gen.is_closure_point()) {
				continue;
			};
			if ( find(gv.begin(), gv.end(), gen) == gv.end() ) {
				gv.push_back(gen);
				Pretropisms.push_back(GeneratorToPoint(gen));
			};
		};
	};
	sort(Pretropisms.begin(), Pretropisms.end());
	PrintPoints(Pretropisms);
	cout << "Number of pretropisms found: " << gv.size() << endl;
	cout << "Hull time: " << HullTime / CLOCKS_PER_SEC << endl;
	cout << "Containment time: " << ContainmentTime / CLOCKS_PER_SEC << endl;
	cout << "Tree append time for loop: " << TreeTime1 / CLOCKS_PER_SEC << endl;
	cout << "Tree append time not in for loop: " << TreeTime2 / CLOCKS_PER_SEC << endl;
	cout << "Intersection time: " << IntersectionTime / CLOCKS_PER_SEC << endl;
	cout << "Parallel time: " << ParallelTime / CLOCKS_PER_SEC << endl;
	cout << "Cleanup time: " << double(clock() - CleanupStart) / CLOCKS_PER_SEC << endl;
	cout << "GetCones time: " << GetConesTime / CLOCKS_PER_SEC << endl;
	if (PreintersectTime != 0) {
		cout << "Preintersection time: " << PreintersectTime / CLOCKS_PER_SEC << endl;
	};
	cout << "Alg time: " << AlgTime / CLOCKS_PER_SEC << endl;
	cout << "Test time: " << TestTime / CLOCKS_PER_SEC << endl;
	cout << "Number of intersections: " << ConeIntersectionCount << endl;
}
