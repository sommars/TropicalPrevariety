#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}


class Edge {
	public:
		C_Polyhedron cone;
		list<int> pointIndices;
		list<int> neighborIndices;
		list<int> childrenIndices;
		list<int> parentIndices;
};

class Facet {
	public:
		list<list<int> > Points;
		C_Polyhedron cone;
		
};

class Hull {
	public:
		list<list<GMP_Integer> > Points;
		map<list<int>,int> PointToIndexMap;
		map<int,list<int> > IndexToPointMap;
		list<Edge> Edges;
		list<Facet> Facets;
		C_Polyhedron CPolyhedron;
};

//------------------------------------------------------------------------------
list<GMP_Integer> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
list<list<GMP_Integer> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
list<GMP_Integer> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(list<list<int> > Points);

//------------------------------------------------------------------------------
int InnerProduct(list<int> V1, list<int> V2);

//------------------------------------------------------------------------------
list<list<int> > FindInitialForm(list<list<int> > Points, list<int> Vector);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(list<list<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(list<int> Point);

//------------------------------------------------------------------------------
void PrintPoint(list<GMP_Integer> Point);
