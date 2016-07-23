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
	Hull NewHull;
	NewHull.CPolyhedron = FindCPolyhedron(Points);
	map<list<GMP_Integer>,GMP_Integer> PointToIndexMap;
	map<GMP_Integer,list<GMP_Integer> > IndexToPointMap;

	NewHull.Points = GeneratorSystemToPoints(NewHull.CPolyhedron.minimized_generators());

	int PtIndex = 0;
	list<list<GMP_Integer> >::iterator itr;
	for (itr=Points.begin(); itr != Points.end(); itr++) {
		list<GMP_Integer>Point=*itr;

		PointToIndexMap[*itr]=PtIndex;
		IndexToPointMap[PtIndex]=*itr;
		PtIndex++;
	}
	NewHull.PointToIndexMap = PointToIndexMap;
	NewHull.IndexToPointMap = IndexToPointMap;
	
	NewHull.Facets = FindFacets(Points, NewHull.CPolyhedron);
	
		//Generators of the lineality space are all of the constraints that are equations.
	
	Constraint_System cs = NewHull.CPolyhedron.minimized_constraints();
	std::cout << "Below are the generators of the lineality space:" << endl;
	for (Constraint_System::const_iterator i = cs.begin(),
		cs_end = cs.end(); i != cs_end; ++i) {

		if ((*i).is_inequality()) {
			std::cout << "Constraints: " << *i << endl;
		}
	}
	
	PrintFacets(NewHull.Facets);
	FindEdges(NewHull);
	return NewHull;
}

//------------------------------------------------------------------------------
list<Facet> FindFacets(list<list<GMP_Integer> > Points, C_Polyhedron ph) {
	// Find the Facets by shooting rays at the polytopes
	// spin through all of the constraints. shoot each constraint at the polytope.
	list<Facet> Facets;
	Constraint_System cs = ph.minimized_constraints();
	cout << "CSSYSTEM" << cs << endl;
	
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		if ((*i).is_inequality()) {
			list<list<GMP_Integer> > FacetPts = FindInitialForm(Points,ConstraintToPoint(*i));	
			Facet NewFacet;
			NewFacet.Points = FacetPts;
			NewFacet.FacetConstraint = *i;
			Facets.push_back(NewFacet);
		};
	}
	return Facets;
}

//------------------------------------------------------------------------------
list<Edge> FindEdges(Hull H) {
//from the set of points, take every pair of points. This is a candidate for being an edge.
//spin through all of the facets. if the candidate edge is sitting on a number of them
//equal to the dimension of the polyhedron - 1, then, the candidate edge is an actual edge.
//the cone of the edge is equal to the C_Polyhedron of all of the constraints of the facets
//as well as the lineality space generators.
	list<Facet> Facets = H.Facets;
	list<list<GMP_Integer> > Points = H.Points;
	list<list<list<GMP_Integer> > > CandidateEdges = FindCandidateEdges(H);
	list<Edge> Edges;
	
	list<list<list<GMP_Integer> > >::iterator itr;
	for (itr=CandidateEdges.begin(); itr != CandidateEdges.end(); itr++) {
		list<list<GMP_Integer> > CandidateEdge = (*itr);
		list<list<GMP_Integer> >::iterator it;
		it=CandidateEdge.begin();
		list<GMP_Integer> Point1 = *it;
		it++;
		list<GMP_Integer> Point2 = *it;
		
		//here it would be good if I could check if Point1 and Point2 are in the set of points on each facet.
		//set<list<
		
		int FacetCount = 0;
		//for facet in facets:
			//	
		
//		if (FacetCount == (possibly polyhedron.dim(). more likely polyhedrom.ambient_dim() - # equations) SUBTRACT 1!!!)
		if (true == true) {
			Edge NewEdge;
			Constraint_System cs;
			Constraint c;
			NewEdge.PointIndices.push_back(H.PointToIndexMap[Point1]);
			NewEdge.PointIndices.push_back(H.PointToIndexMap[Point2]);
			NewEdge.Cone = C_Polyhedron(cs);
			Edges.push_back(NewEdge);
		}
	};
	
	
	// After all of the edges have been generated, fill out all of the neighbors on all of the edges.

	
	return Edges;
}

//------------------------------------------------------------------------------
list<list<list<GMP_Integer> > > FindCandidateEdges(Hull H) {


	list<list<list<GMP_Integer> > > CandidateEdges;
	int n = H.Points.size();
	vector<int> d(n);
	for (size_t i = 0; i != d.size(); ++i) {
		d[i] = i;
	}
	do {
		list<list<GMP_Integer> > CandidateEdge;
		for (int i = 0; i < 2; i++) {
			CandidateEdge.push_back(H.IndexToPointMap[d[i]]);
		}
		CandidateEdges.push_back(CandidateEdge);
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
