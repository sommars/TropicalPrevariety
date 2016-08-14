#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>
#include <string>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron &ph1, C_Polyhedron &ph2) {
	Constraint_System cs1 = ph1.minimized_constraints();
	Constraint_System cs2 = ph2.minimized_constraints();
	for (Constraint_System::const_iterator i = cs1.begin(),
	cs1_end = cs1.end(); i != cs1_end; ++i) {
		cs2.insert(*i);
	}
	Recycle_Input dummy;
	C_Polyhedron ph(cs2,dummy);
	ph.affine_dimension();
	return ph;
}

//------------------------------------------------------------------------------
Cone IntersectCones(Cone C1, Cone C2) {
	Constraint_System cs1 = C1.Polyhedron.minimized_constraints();
	Constraint_System cs2 = C2.Polyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs1.begin(),
	cs1_end = cs1.end(); i != cs1_end; ++i) {
		cs2.insert(*i);
	}
	Recycle_Input dummy;
	C_Polyhedron ph(cs2,dummy);
	ph.affine_dimension();
	Cone C3;
	C3.Polyhedron = ph;
	return C3;
}

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < g.space_dimension(); i++) {
		stringstream s;
		s << (g).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
vector<int> GeneratorToIntPoint(Generator g) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < g.space_dimension(); i++) {
		stringstream s;
		s << (g).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
vector<vector<int> > GeneratorSystemToPoints(Generator_System gs) {
	vector<vector<int> > Result;
	for (Generator_System::const_iterator i = gs.begin(),
	gs_end = gs.end(); i != gs_end; ++i) {
		Result.push_back(GeneratorToPoint(*i));
	}
	return Result;
}

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c) { //page 251
	vector<int> Result;
	for (size_t i = 0; i < c.space_dimension(); i++) {
		stringstream s;
		s << (c).coefficient(Variable(i));
		int ToAppend;
		istringstream(s.str()) >> ToAppend;
		Result.push_back(ToAppend);
	}

	return Result;
}

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<int> > Points) {
	Hull H;

	H.CPolyhedron = FindCPolyhedron(Points);
	H.Points = GeneratorSystemToPoints(H.CPolyhedron.minimized_generators());
	H.AffineDimension = H.CPolyhedron.affine_dimension();
	H.SpaceDimension = H.CPolyhedron.space_dimension();
	
	// Find the lineality space.
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		Constraint c = *i;
		if (c.is_equality()) {
			Linear_Expression LE;
			for (dimension_type iii = c.space_dimension(); iii-- > 0; ) {
				LE += c.coefficient(Variable(iii)) * Variable(iii);
			};
			H.LinealitySpace.insert(line(LE));
		};
	};
	
	// Create PointToIndexMap
	for (size_t i = 0; i != H.Points.size(); i++) {
		vector<int> Point = H.Points[i];
		H.PointToIndexMap[Point]=i;
	};

	H.Facets = FindFacets(H);
	H.Edges = FindEdges(H);
	cout << "Convex hull------------------------" << endl;
	PrintPoints(H.Points);
	cout << "Affine dimension: " << H.AffineDimension << endl;
	cout << "Space dimension: " << H.SpaceDimension << endl;
	cout << "Number of edges: " << H.Edges.size() << endl;
	cout << "Number of facets: " << H.Facets.size() << endl << endl;
	return H;
}

//------------------------------------------------------------------------------
vector<Facet> FindFacets(Hull H) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	C_Polyhedron ph = H.CPolyhedron;
	vector<vector<int> > Points = H.Points;
	vector<Facet> Facets;
	Constraint_System cs = ph.minimized_constraints();
	
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		if (!(*i).is_inequality()) {
			continue;
		};
		vector<int> Pt = ConstraintToPoint(*i);
		vector<vector<int> > FacetPts = FindInitialForm(Points, Pt);
		Facet F;
		Generator_System gs;
		vector<int>::iterator it;
		Linear_Expression RayLE;
		Linear_Expression PointLE;
		for(size_t VarIndex = 0; VarIndex != Pt.size(); VarIndex++) {
			RayLE += Variable(VarIndex) * (Pt[VarIndex]);
			PointLE += Variable(VarIndex) * 0;
		};
		gs.insert(ray(RayLE));
		gs.insert(point(PointLE));

		F.Generators = C_Polyhedron(gs).minimized_generators();
		
		vector<vector<int> >::iterator itr;
		
		for (itr=FacetPts.begin(); itr != FacetPts.end(); itr++) {
			F.PointIndices.insert(H.PointToIndexMap[*itr]);
		};
		Facets.push_back(F);
	}
	return Facets;
}

//------------------------------------------------------------------------------
vector<Edge> FindEdges(Hull H) {
	vector<Facet> Facets = H.Facets;
	vector<vector<int> > Points = H.Points;
	vector<vector<int> > CandidateEdges = FindCandidateEdges(H);
	vector<Edge> Edges;

	//The number of facets we want is equal to the dimension of the ambient space minus the number of equations -1
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	int Dim = H.CPolyhedron.affine_dimension() - 1;

	vector<vector<int> >::iterator itr;
	for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++) {
		vector<int> CandidateEdge = (*itr);
		int Point1 = CandidateEdge[0];
		int Point2 = CandidateEdge[1];

		int FacetCount = 0;
		Generator_System gs;
		for (vector<Facet>::iterator FacetIt=Facets.begin(); FacetIt != Facets.end(); FacetIt++) {
			set<int> PtIndices = (*FacetIt).PointIndices;
			bool Point1IsInFacet = PtIndices.find(Point1) != PtIndices.end();
			bool Point2IsInFacet = PtIndices.find(Point2) != PtIndices.end();
			if (Point1IsInFacet and Point2IsInFacet) {
				Generator_System Generators = (*FacetIt).Generators;
				for (Generator_System::const_iterator GenItr = Generators.begin(),
				cs_end = Generators.end(); GenItr != cs_end; ++GenItr) {
					gs.insert(*GenItr);
				};
				FacetCount++;
			};
		};

		if (FacetCount >= Dim) {
			Edge NewEdge;
			NewEdge.PointIndices.insert(Point1);
			NewEdge.PointIndices.insert(Point2);

			// Need to add all of the generators of the lineality space to the cones.
			Generator_System gs2 = H.LinealitySpace;
			for (Generator_System::const_iterator i = gs2.begin(),
			cs_end = gs2.end(); i != cs_end; ++i) {
				gs.insert(*i);
			};
			NewEdge.EdgeCone.Polyhedron = C_Polyhedron(gs);
			
			// Force PPL to calculate as much as we need ahead of time.
			NewEdge.EdgeCone.Polyhedron.minimized_generators();
			NewEdge.EdgeCone.Polyhedron.minimized_constraints();
			NewEdge.EdgeCone.Polyhedron.affine_dimension();
			Edges.push_back(NewEdge);
		}
	};

	// After all of the edges have been generated, fill out all of the neighbors on all of the edges.
	for (size_t Edge1Index = 0; Edge1Index != Edges.size(); Edge1Index++) {
		Edge Edge1 = Edges[Edge1Index];

		for (size_t Edge2Index = 0; Edge2Index != Edges.size(); Edge2Index++) {
			if (Edge1Index == Edge2Index) {
				continue;
			};

			int IntersectionCount = IntersectSets(Edge1.PointIndices, Edges[Edge2Index].PointIndices).size();
			if (IntersectionCount == 1) {
				Edges[Edge1Index].NeighborIndices.insert(Edge2Index);
				Edges[Edge2Index].NeighborIndices.insert(Edge1Index);
			} else if (IntersectionCount > 1) {
				cout << "Internal Error: Two edges intersect at more than one point" << endl;
				cin.get();
			}
		};
	};

	return Edges;
}

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull H) {
	vector<vector<int> > CandidateEdges;
	int n = H.Points.size();
	vector<int> d(n);
	for (size_t i = 0; i != d.size(); ++i) {
		d[i] = i;
	}
	do {
		if (d[0] < d[1]) {
			vector<int> CandidateEdgeIndices;
			CandidateEdgeIndices.push_back(d[0]);
			CandidateEdgeIndices.push_back(d[1]);
			CandidateEdges.push_back(CandidateEdgeIndices);
		};
		reverse(d.begin()+2, d.end());
	} while (next_permutation(d.begin(), d.end()));
	return CandidateEdges;
}

//------------------------------------------------------------------------------
int InnerProduct(vector<int> V1, vector<int> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
		cin.get();
	};
	int Result = 0;
	for (size_t i = 0; i != V1.size(); i++) {
		Result += V1[i] * V2[i];
	}
	return Result;
}

//------------------------------------------------------------------------------
vector<vector<int> > FindInitialForm(vector<vector<int> > &Points, vector<int> &Vector) {
	/*
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};
	vector<vector<int> > InitialForm;

	InitialForm.push_back(Points[0]);
	int MinimalIP = InnerProduct(Vector, Points[0]);

	for (size_t i = 1; i != Points.size(); i++) {
		vector<int> Point = Points[i];
		int IP = InnerProduct(Vector, Point);
		if (MinimalIP > IP) {
			MinimalIP = IP;
			InitialForm.clear();
			InitialForm.push_back(Point);
		} else if (IP == MinimalIP) {
			InitialForm.push_back(Point);
		};
	};

	return InitialForm;
}

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points) {
	Generator_System gs;
	vector<vector<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		vector<int> Point = *itr;
		Linear_Expression LE;
		for (size_t i = 0; i != Point.size(); i++) {
			LE = LE + Variable(i) * (Point[i]);
		};
		gs.insert(point(LE));
	};
	C_Polyhedron ph = C_Polyhedron(gs);
	
	return ph;
}

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point) {
	vector<int>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points) {
	vector<vector<int> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		PrintPoint(*itr);
	};
}

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point) {
	set<int>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintCPolyhedron(C_Polyhedron ph, bool PrintIf0Dim) {
	cout << "C_Polyhedron dimension:" << ph.affine_dimension() << endl;
	cout << "Generators: " << ph.minimized_generators() << endl;
}

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<C_Polyhedron> phs, bool PrintIf0Dim) {
	vector<C_Polyhedron>::iterator it;
	for (it=phs.begin(); it != phs.end(); it++) {
		PrintCPolyhedron(*it, PrintIf0Dim);
	};
}

//------------------------------------------------------------------------------
set<int> IntersectSets(set<int> S1, set<int> S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	set<int> Result;
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			Result.insert(*S1Itr);
			S1Itr++;
			S2Itr++;
		};
	};
	return Result;
}

//------------------------------------------------------------------------------
bool SetsDoIntersect(set<int> &S1, set<int> &S2) {
	set<int>::iterator S1Itr = S1.begin();
	set<int>::iterator S2Itr = S2.begin();
	while ((S1Itr != S1.end()) && (S2Itr != S2.end())) {
		if (*S1Itr < *S2Itr) {
			++S1Itr;
		}
		else if (*S2Itr<*S1Itr) {
			++S2Itr;
		} else {
			return true;
		};
	};
	return false;
}
