#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}
double IntersectionTime;

//------------------------------------------------------------------------------
double GetPolyhedralIntersectionTime() {
	return IntersectionTime;
}

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron ph1, C_Polyhedron ph2) {
	
  clock_t begin = clock();
	Constraint_System cs;
	
	Constraint_System cs1 = ph1.minimized_constraints();
	Constraint_System cs2 = ph2.minimized_constraints();
	for (Constraint_System::const_iterator i = cs1.begin(),
	cs1_end = cs1.end(); i != cs1_end; ++i) {
		cs.insert(*i);
	}
	for (Constraint_System::const_iterator i = cs2.begin(),
	cs2_end = cs2.end(); i != cs2_end; ++i) {
		cs.insert(*i);
	}
	C_Polyhedron ph(cs);
	ph.minimized_constraints();

	IntersectionTime += double(clock() - begin);
	return ph;
}

//------------------------------------------------------------------------------
vector<GMP_Integer> GeneratorToPoint(Generator g) { //page 251
	vector<GMP_Integer> Result;
	for (int i = 0; i < g.space_dimension(); i++) {
		Result.push_back((g).coefficient(Variable(i)));
	}

	return Result;
}

//------------------------------------------------------------------------------
vector<GMP_Integer> ParseFinalOutput(Generator g) { //page 251
	vector<GMP_Integer> Result;
	bool reverse = false;
	
	for (int i = 0; i < g.space_dimension(); i++) {
		if (g.coefficient(Variable(i)) != 0) {
			if (g.coefficient(Variable(i)) < 0){
				reverse = true;
			};
			break;
		};
	}	
	for (int i = 0; i < g.space_dimension(); i++) {
		Result.push_back((g).coefficient(Variable(i)));
	}

	return Result;
}

//------------------------------------------------------------------------------
vector<vector<GMP_Integer> > GeneratorSystemToPoints(Generator_System gs) {
	vector<vector<GMP_Integer> > Result;

	for (Generator_System::const_iterator i = gs.begin(),
	gs_end = gs.end(); i != gs_end; ++i) {
		Result.push_back(GeneratorToPoint(*i));
	}
	return Result;
}

//------------------------------------------------------------------------------
vector<GMP_Integer> ConstraintToPoint(Constraint c) { //page 251
	vector<GMP_Integer> Result;
	for (int i = 0; i < c.space_dimension(); i++) {
		Result.push_back((c).coefficient(Variable(i)));
	}
	//Note: This is commented out for the sake of FindFacets
	//Result.push_back(c.inhomogeneous_term());
	return Result;
}

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<GMP_Integer> > Points) {
	Hull H;

	H.CPolyhedron = FindCPolyhedron(Points);
	H.Points = GeneratorSystemToPoints(H.CPolyhedron.minimized_generators());
	H.Dimension = H.CPolyhedron.affine_dimension();
	
	// Find the lineality space.
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	Generator_System gs;
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		Constraint c = *i;
		if (c.is_equality()) {
			Linear_Expression LE;
			for (dimension_type iii = c.space_dimension(); iii-- > 0; ) {
				LE += c.coefficient(Variable(iii)) * Variable(iii);
			};
			gs.insert(line(LE));
		};
	};
	H.LinealitySpace = gs;
	
	// Create point/index maps.
	map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap;
	map<GMP_Integer,vector<GMP_Integer> > IndexToPointMap;
	for (int i = 0; i != H.Points.size(); i++) {
		vector<GMP_Integer> Point = H.Points[i];
		PointToIndexMap[Point]=i;
		IndexToPointMap[i]=Point;
	};
	H.PointToIndexMap = PointToIndexMap;
	H.IndexToPointMap = IndexToPointMap;

	H.Facets = FindFacets(H);
	H.Edges = FindEdges(H);
	cout << "Convex hull------------------------" << endl;
	PrintPoints(H.Points);
	cout << "Dimension: " << H.Dimension << endl;
	cout << "Number of edges: " << H.Edges.size() << endl;
	cout << "Number of facets: " << H.Facets.size() << endl << endl;
	return H;
}

//------------------------------------------------------------------------------
vector<Facet> FindFacets(Hull H) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	C_Polyhedron ph = H.CPolyhedron;
	vector<vector<GMP_Integer> > Points = H.Points;
	vector<Facet> Facets;
	Constraint_System cs = ph.minimized_constraints();
	
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		if ((*i).is_inequality()) {
			vector<GMP_Integer> Pt = ConstraintToPoint(*i);
			vector<vector<GMP_Integer> > FacetPts = FindInitialForm(Points, Pt);
			Facet F;
			F.Points = FacetPts;
			F.FacetConstraint = *i;
			Generator_System gs;
			vector<GMP_Integer>::iterator it;
			Linear_Expression LE1;
			Linear_Expression LE2;
			int VarIndex = 0;
			for (it=Pt.begin(); it != Pt.end(); it++) {
				LE1 = LE1 + Variable(VarIndex) * (*it);
				LE2 = LE2 + Variable(VarIndex) * 0;
				VarIndex++;
			};
			gs.insert(ray(LE1));
			gs.insert(point(LE2));
			C_Polyhedron CPoly = C_Polyhedron(gs);
			
			F.ConeConstraints = CPoly.minimized_constraints();
			
			set<GMP_Integer> PointIndices;
			
			vector<vector<GMP_Integer> >::iterator itr;
			
			for (itr=FacetPts.begin(); itr != FacetPts.end(); itr++) {
				PointIndices.insert(H.PointToIndexMap[*itr]);
			};
			F.PointIndices = PointIndices;
			Facets.push_back(F);
		};
	}
	return Facets;
}

//------------------------------------------------------------------------------
vector<Edge> FindEdges(Hull H) {
	vector<Facet> Facets = H.Facets;
	vector<vector<GMP_Integer> > Points = H.Points;
	vector<vector<GMP_Integer> > CandidateEdges = FindCandidateEdges(H);
	vector<Edge> Edges;

	//The number of facets we want is equal to the dimension of the ambient space minus the number of equations -1
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	int Dim = H.CPolyhedron.affine_dimension() - 1;

	vector<vector<GMP_Integer> >::iterator itr;
	for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++) {
		vector<GMP_Integer> CandidateEdge = (*itr);
		vector<GMP_Integer>::iterator it;
		it=CandidateEdge.begin();
		GMP_Integer Point1 = *it;
		it++;
		GMP_Integer Point2 = *it;

		int FacetCount = 0;
		vector<Facet>::iterator FacetIt;

		Generator_System gs;
		for (FacetIt=Facets.begin(); FacetIt != Facets.end(); FacetIt++) {
			set<GMP_Integer> PtIndices = (*FacetIt).PointIndices;
			bool Point1IsInFacet = PtIndices.find(Point1) != PtIndices.end();
			bool Point2IsInFacet = PtIndices.find(Point2) != PtIndices.end();
			if (Point1IsInFacet and Point2IsInFacet) {
				FacetCount++;
				Generator_System ConeGenerators = C_Polyhedron((*FacetIt).ConeConstraints).minimized_generators();
				for (Generator_System::const_iterator iii = ConeGenerators.begin(),
				cs_end = ConeGenerators.end(); iii != cs_end; ++iii) {
					gs.insert(*iii);
				}
			};
		};

		// Is >= correct?
		if (FacetCount >= Dim) {
			Edge NewEdge;
			set<GMP_Integer> PointIndices;
			set<int> NeighborIndices;
			PointIndices.insert(Point1);
			PointIndices.insert(Point2);
			NewEdge.PointIndices = PointIndices;
			NewEdge.NeighborIndices = NeighborIndices;

			// Need to add all of the generators of the lineality space to the cones.
			Generator_System gs2 = H.LinealitySpace;
			for (Generator_System::const_iterator i = gs2.begin(),
			cs_end = gs2.end(); i != cs_end; ++i) {
				gs.insert(*i);
			};
			NewEdge.Cone = C_Polyhedron(gs);
			Edges.push_back(NewEdge);
		}
	};
	
	// After all of the edges have been generated, fill out all of the neighbors on all of the edges.
	for (int Edge1Index = 0; Edge1Index != Edges.size(); Edge1Index++) {
		Edge Edge1 = Edges[Edge1Index];
		set<GMP_Integer>::iterator SetItr1=Edge1.PointIndices.begin();
		GMP_Integer PtIndex1 = *SetItr1;
		SetItr1++;
		GMP_Integer PtIndex2 = *SetItr1;
		
		for (int Edge2Index = 0; Edge2Index != Edges.size(); Edge2Index++) {
			if (Edge1Index == Edge2Index) {
				continue;
			};
			Edge Edge2 = Edges[Edge2Index];
			
			int IntersectionCount = 0;
			for (std::set<GMP_Integer>::iterator SetItr2=Edge2.PointIndices.begin(); 
			SetItr2!=Edge2.PointIndices.end(); ++SetItr2) {
				if ((*SetItr2 == PtIndex1) or (*SetItr2 == PtIndex2)) {
					IntersectionCount++;
				}
			}
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
vector<vector<GMP_Integer> > FindCandidateEdges(Hull H) {
	vector<vector<GMP_Integer> > CandidateEdges;
	int n = H.Points.size();
	vector<int> d(n);
	for (size_t i = 0; i != d.size(); ++i) {
		d[i] = i;
	}
	do {
		if (d[0] < d[1]) {
			vector<GMP_Integer> CandidateEdgeIndices;
			CandidateEdgeIndices.push_back(d[0]);
			CandidateEdgeIndices.push_back(d[1]);
			CandidateEdges.push_back(CandidateEdgeIndices);
		};
		reverse(d.begin()+2, d.end());
	} while (next_permutation(d.begin(), d.end()));
	return CandidateEdges;
}

//------------------------------------------------------------------------------
GMP_Integer InnerProduct(vector<GMP_Integer> V1, vector<GMP_Integer> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
		cin.get();
	};
	GMP_Integer Result = 0;
	for (int i = 0; i != V1.size(); i++) {
		Result += V1[i] * V2[i];
	}
	return Result;
}

//------------------------------------------------------------------------------
vector<vector<GMP_Integer> > FindInitialForm(vector<vector<GMP_Integer> > Points, vector<GMP_Integer> Vector) {
	/*
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};
	vector<vector<GMP_Integer> > InitialForm;

	InitialForm.push_back(Points[0]);
	GMP_Integer MinimalIP = InnerProduct(Vector, Points[0]);

	for (int i = 1; i != Points.size(); i++) {
		vector<GMP_Integer> Point = Points[i];
		GMP_Integer IP = InnerProduct(Vector, Point);
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
C_Polyhedron FindCPolyhedron(vector<vector<GMP_Integer> > Points) {
	Generator_System gs;
	vector<vector<GMP_Integer> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		vector<GMP_Integer> Point = *itr;
		Linear_Expression LE;
		for (int i = 0; i != Point.size(); i++) {
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
void PrintPoint(vector<GMP_Integer> Point) {
	vector<GMP_Integer>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<GMP_Integer> > Points) {
	vector<vector<GMP_Integer> >::iterator itr;
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
set<GMP_Integer> IntersectSets(set<GMP_Integer> S1, set<GMP_Integer> S2) {
	set<GMP_Integer> Result;
	set<GMP_Integer>::iterator SetItr;
	for (SetItr = S1.begin(); SetItr != S1.end(); SetItr++) {
		GMP_Integer TestInt = *SetItr;
		if ( find(S2.begin(), S2.end(), TestInt) != S2.end() ) {
			if ( find(Result.begin(), Result.end(), TestInt) == Result.end() ) {
				Result.insert(TestInt);
			};
		};
	};
	return Result;
}
