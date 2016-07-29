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
		set<int> NeighborIndices;
};

class Facet {
	public:
		set<GMP_Integer> PointIndices;
		Generator_System Generators;
};

class Hull {
	public:
		vector<vector<GMP_Integer> > Points;
		map<vector<GMP_Integer>,GMP_Integer> PointToIndexMap;
		vector<Edge> Edges;
		vector<Facet> Facets;
		C_Polyhedron CPolyhedron;
		Generator_System LinealitySpace;
		int Dimension;
};

//------------------------------------------------------------------------------
double GetPolyhedralIntersectionTime();

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron ph1, C_Polyhedron ph2);

//------------------------------------------------------------------------------
vector<GMP_Integer> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
vector<vector<GMP_Integer> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
vector<GMP_Integer> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<GMP_Integer> > Points);

//------------------------------------------------------------------------------
vector<Facet> FindFacets(Hull H);

//------------------------------------------------------------------------------
vector<Edge> FindEdges(Hull H);

//------------------------------------------------------------------------------
vector<vector<GMP_Integer> > FindCandidateEdges(Hull H);

//------------------------------------------------------------------------------
GMP_Integer InnerProduct(vector<GMP_Integer> V1, vector<GMP_Integer> V2);

//------------------------------------------------------------------------------
vector<vector<GMP_Integer> > FindInitialForm(vector<vector<GMP_Integer> > Points, vector<GMP_Integer> Vector);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<GMP_Integer> > Points);

//------------------------------------------------------------------------------
void PrintPoint(vector<int> Point);

//------------------------------------------------------------------------------
void PrintPoint(vector<GMP_Integer> Point);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<GMP_Integer> > Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point);

//------------------------------------------------------------------------------
void PrintCPolyhedron(C_Polyhedron ph, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<C_Polyhedron> phs, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
set<GMP_Integer> IntersectSets(set<GMP_Integer> S1, set<GMP_Integer> S2);

//------------------------------------------------------------------------------
bool SetsDoIntersect(set<GMP_Integer> S1, set<GMP_Integer> S2);
