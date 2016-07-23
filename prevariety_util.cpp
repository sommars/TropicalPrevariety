#include "prevariety_util.h"
#include <iostream>
#include <list>
#include <stdio.h>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//------------------------------------------------------------------------------
list<GMP_Integer> GeneratorToPoint(Generator g) { //page 251
	list<GMP_Integer> Result;
	for (int i = 0; i < g.space_dimension(); i++) {
		Result.push_back((g).coefficient(Variable(i)));
	}

	return Result;
}

//------------------------------------------------------------------------------
list<list<GMP_Integer> > GeneratorSystemToPoints(Generator_System gs) {
	list<list<GMP_Integer> > Result;

	for (Generator_System::const_iterator i = gs.begin(),
	gs_end = gs.end(); i != gs_end; ++i) {
		Result.push_back(GeneratorToPoint(*i));
	}
	return Result;
}

//------------------------------------------------------------------------------
list<GMP_Integer> ConstraintToPoint(Constraint c) { //page 251
	list<GMP_Integer> Result;
	for (int i = 0; i < c.space_dimension(); i++) {
		Result.push_back((c).coefficient(Variable(i)));
	}
	//Note: This is commented out for the sake of FindFacets
	//Result.push_back(c.inhomogeneous_term());
	return Result;
}

//------------------------------------------------------------------------------
Hull NewHull(list<list<GMP_Integer> > Points) {
	Hull H;
	// Find the C_Polyhedron
	H.CPolyhedron = FindCPolyhedron(Points);

	H.Points = GeneratorSystemToPoints(H.CPolyhedron.minimized_generators());

	// Create point/index maps.
	map<list<GMP_Integer>,GMP_Integer> PointToIndexMap;
	map<GMP_Integer,list<GMP_Integer> > IndexToPointMap;
	int PtIndex = 0;
	list<list<GMP_Integer> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		list<GMP_Integer>Point=*itr;

		PointToIndexMap[*itr]=PtIndex;
		IndexToPointMap[PtIndex]=*itr;
		PtIndex++;
	}
	H.PointToIndexMap = PointToIndexMap;
	H.IndexToPointMap = IndexToPointMap;
	
	// Find the facets.
	H.Facets = FindFacets(H);
	
	// Find the lineality space.
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	Constraint_System LS;
	std::cout << "Below are the generators of the lineality space:" << endl;
	for (Constraint_System::const_iterator i = cs.begin(),
		cs_end = cs.end(); i != cs_end; ++i) {
		if ((*i).is_inequality()) {
			LS.insert(*i);
		};
	};
	H.LinealitySpace = LS;
	
	// Find the edges.
	H.Edges = FindEdges(H);
	return H;
}

//------------------------------------------------------------------------------
list<Facet> FindFacets(Hull H) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	C_Polyhedron ph = H.CPolyhedron;
	list<list<GMP_Integer> > Points = H.Points;
	list<Facet> Facets;
	Constraint_System cs = ph.minimized_constraints();
	cout << "CSSYSTEM" << cs << endl;
	
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		if ((*i).is_inequality()) {
			list<list<GMP_Integer> > FacetPts = FindInitialForm(Points,ConstraintToPoint(*i));	
			Facet F;
			F.Points = FacetPts;
			F.FacetConstraint = *i;
			set<GMP_Integer> PointIndices;
			
			list<list<GMP_Integer> >::iterator itr;
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
list<Edge> FindEdges(Hull H) {
	list<Facet> Facets = H.Facets;
	list<list<GMP_Integer> > Points = H.Points;
	list<list<GMP_Integer> > CandidateEdges = FindCandidateEdges(H);
	list<Edge> Edges;

	//The number of facets we want is equal to the dimension of the ambient space minus the number of equations -1
	Constraint_System cs = H.CPolyhedron.minimized_constraints();
	int EquationCount = 0;	
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		if ((*i).is_equality()) {
			EquationCount++;
		};
	};
	int Dim = H.CPolyhedron.space_dimension() - EquationCount - 1;

	list<list<GMP_Integer> >::iterator itr;
	for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++) {
		list<GMP_Integer> CandidateEdge = (*itr);
		list<GMP_Integer>::iterator it;
		it=CandidateEdge.begin();
		GMP_Integer Point1 = *it;
		it++;
		GMP_Integer Point2 = *it;

		int FacetCount = 0;
		list<Facet>::iterator FacetIt;
		Constraint_System cs;
		for (FacetIt=Facets.begin(); FacetIt != Facets.end(); it++) {
			set<GMP_Integer> PtIndices = (*FacetIt).PointIndices;
			bool Point1IsInFacet = PtIndices.find(Point1) != PtIndices.end();
			bool Point2IsInFacet = PtIndices.find(Point2) != PtIndices.end();
			if (Point1IsInFacet and Point2IsInFacet) {
				FacetCount++;
				cs.insert((*FacetIt).FacetConstraint);
			};
		};

		if (FacetCount == Dim) {
			Edge NewEdge;
			set<GMP_Integer> PointIndices;
			set<GMP_Integer> NeighborIndices;
			PointIndices.insert(Point1);
			PointIndices.insert(Point2);
			NewEdge.PointIndices = PointIndices;
			NewEdge.NeighborIndices = NeighborIndices;

			// Need to add all of the generators of the lineality space to the cones.
			for (Constraint_System::const_iterator i = H.LinealitySpace.begin(),
			cs_end = H.LinealitySpace.end(); i != cs_end; ++i) {
				cs.insert(*i);
			};
			NewEdge.Cone = C_Polyhedron(cs);
			Edges.push_back(NewEdge);
		}
	};
	
	// After all of the edges have been generated, fill out all of the neighbors on all of the edges.
	
	int Edge1Index = 0;

	list<Edge>::iterator EdgeItr1;
	for (EdgeItr1=Edges.begin(); EdgeItr1 != Edges.end(); EdgeItr1++) {
		list<Edge>::iterator EdgeItr2;
		int Edge2Index = 0;
		Edge Edge1 = *EdgeItr1;
		
		set<GMP_Integer>::iterator SetItr1=Edge1.PointIndices.begin();
		GMP_Integer PtIndex1 = *SetItr1;
		SetItr1++;
		GMP_Integer PtIndex2 = *SetItr1;

		for (EdgeItr2=Edges.begin(); EdgeItr2 != Edges.end(); EdgeItr2++) {
			if (Edge1Index == Edge2Index) {
				Edge2Index++;
				continue;
			};
			Edge Edge2 = *EdgeItr2;

			int IntersectionCount = 0;
			for (std::set<GMP_Integer>::iterator SetItr2=Edge2.PointIndices.begin(); 
			SetItr2!=Edge2.PointIndices.end(); ++SetItr2) {
				if ((*SetItr2 == PtIndex1) or (*SetItr2 == PtIndex2)) {
					IntersectionCount++;
				}
			}

			//if the intersection is one, then the neighbors are edges.
			if (IntersectionCount == 1) {
				Edge1.NeighborIndices.insert(Edge2Index);
				Edge2.NeighborIndices.insert(Edge1Index);
			} else if (IntersectionCount > 1) {
				cout << "Internal Error: Two edges intersect at more than one point" << endl;
			}
			Edge2Index++;
		}
		Edge1Index++;	
	};
	return Edges;
}

//------------------------------------------------------------------------------
list<list<GMP_Integer> > FindCandidateEdges(Hull H) {

	list<list<GMP_Integer> > CandidateEdges;
	int n = H.Points.size();
	vector<int> d(n);
	for (size_t i = 0; i != d.size(); ++i) {
		d[i] = i;
	}
	do {
		list<GMP_Integer> CandidateEdgeIndices;
		for (int i = 0; i < 2; i++) {
			CandidateEdgeIndices.push_back(d[i]);
		}
		CandidateEdges.push_back(CandidateEdgeIndices);
		reverse(d.begin()+2, d.end());
	} while (next_permutation(d.begin(), d.end()));
	return CandidateEdges;
}

//------------------------------------------------------------------------------
GMP_Integer InnerProduct(list<GMP_Integer> V1, list<GMP_Integer> V2) {
	/* 
		Computes the inner product of two vectors.
	*/
	if (V1.size() != V2.size()) {
		cout << "Internal Error: InnerProduct with different sizes" << endl;
	};
	GMP_Integer Result = 0;
	list<GMP_Integer>::iterator it1;
	list<GMP_Integer>::iterator it2;
	it2 = V2.begin();
	for (it1=V1.begin(); it1 != V1.end(); it1++) {
		Result = Result + (*it1) * (*it2);
		it2++;
	}
	return Result;
}

//------------------------------------------------------------------------------
list<list<GMP_Integer> > FindInitialForm(list<list<GMP_Integer> > Points, list<GMP_Integer> Vector) {
	/* 
		Computes the initial form of a vector and a set of points.
	*/
	if (Points.size() == 0) {
		return Points;
	};
	list<list<GMP_Integer> > IF;
	list<list<GMP_Integer> >::iterator itr;

	itr=Points.begin();
	IF.push_back(*itr);
	GMP_Integer MinimalIP = InnerProduct(Vector, *itr);
	itr++;

	for (itr; itr != Points.end(); itr++) {
		GMP_Integer IP = InnerProduct(Vector, *itr);
		if (MinimalIP > IP) {
			MinimalIP = IP;
			IF.clear();
			IF.push_back(*itr);
		} else if (IP == MinimalIP) {
			IF.push_back(*itr);
		}
	}
	
	return IF;
}

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(list<list<GMP_Integer> > Points) {
	Generator_System gs;
	list<list<GMP_Integer> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		list<GMP_Integer>Point=*itr;
		list<GMP_Integer>::iterator it;
		Linear_Expression LE;
		int VarIndex = 0;
		for (it=Point.begin(); it != Point.end(); it++) {
			LE = LE + Variable(VarIndex) * (*it);
			VarIndex++;
		}
		gs.insert(point(LE));
	}
	C_Polyhedron ph = C_Polyhedron(gs);
	return ph;
}

//------------------------------------------------------------------------------
void PrintPoint(list<int> Point) {
	list<int>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoint(list<GMP_Integer> Point) {
	list<GMP_Integer>::iterator it;
	cout << "{ ";
	for (it=Point.begin(); it != Point.end(); it++) {
		cout << (*it) << " ";
	}
	cout << "}" << endl;
}

//------------------------------------------------------------------------------
void PrintPoints(list<list<GMP_Integer> > Points) {
	list<list<GMP_Integer> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		PrintPoint(*itr);
	};
}

//------------------------------------------------------------------------------
void PrintFacet(Facet F){
	cout << "FacetPts below"<< endl;
	PrintPoints(F.Points);
	cout << endl;
	cout << "Constraint: " << F.FacetConstraint << endl;
}

//------------------------------------------------------------------------------
void PrintFacets(list<Facet> Facets) {
	list<Facet>::iterator it;
	cout << "Printing Facets-------------------------------" << endl;
	for (it=Facets.begin(); it != Facets.end(); it ++) {
		cout << "NEW FACET" << endl;
	};
}
