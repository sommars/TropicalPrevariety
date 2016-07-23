#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}


class Edge {
	public:
		C_Polyhedron Cone;
		set<GMP_Integer> PointIndices;
		set<GMP_Integer> NeighborIndices;
};

class Facet {
	public:
		list<list<GMP_Integer> > Points; //this should be a set
		set<GMP_Integer> PointIndices;
		Constraint FacetConstraint;
		Constraint_System ConeConstraints;
};

class Hull {
	public:
		list<list<GMP_Integer> > Points; //this should be a set
		map<list<GMP_Integer>,GMP_Integer> PointToIndexMap;
		map<GMP_Integer,list<GMP_Integer> > IndexToPointMap;
		list<Edge> Edges;
		list<Facet> Facets;
		C_Polyhedron CPolyhedron;
		Constraint_System LinealitySpace;
};

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron ph1, C_Polyhedron ph2);

//------------------------------------------------------------------------------
list<GMP_Integer> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
list<list<GMP_Integer> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
list<GMP_Integer> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(list<list<GMP_Integer> > Points);

//------------------------------------------------------------------------------
list<Facet> FindFacets(Hull H);

//------------------------------------------------------------------------------
list<Edge> FindEdges(Hull H);

//------------------------------------------------------------------------------
list<list<GMP_Integer> > FindCandidateEdges(Hull H);

//------------------------------------------------------------------------------
GMP_Integer InnerProduct(list<GMP_Integer> V1, list<GMP_Integer> V2);

//------------------------------------------------------------------------------
list<list<GMP_Integer> > FindInitialForm(list<list<GMP_Integer> > Points, list<GMP_Integer> Vector);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(list<list<GMP_Integer> > Points);

//------------------------------------------------------------------------------
void PrintPoint(list<int> Point);

//------------------------------------------------------------------------------
void PrintPoint(list<GMP_Integer> Point);

//------------------------------------------------------------------------------
void PrintPoints(list<list<GMP_Integer> > Points);

//------------------------------------------------------------------------------
void PrintFacet(Facet F);

//------------------------------------------------------------------------------
void PrintFacets(list<Facet> Facets);

//------------------------------------------------------------------------------
int FindCSDim(Constraint_System cs);

//------------------------------------------------------------------------------
void PrintCPolyhedron(C_Polyhedron ph, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
void PrintCPolyhedrons(list<C_Polyhedron> phs, bool PrintIf0Dim = true);
