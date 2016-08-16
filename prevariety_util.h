#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

class Cone {
	public:
		C_Polyhedron Polyhedron;
		vector<set<int> > IntersectionIndices;
};

class Edge {
	public:
		Cone EdgeCone;
		set<int> PointIndices;
		set<int> NeighborIndices;
};

class Facet {
	public:
		set<int> PointIndices;
		Generator_System Generators;
};

class Hull {
	public:
		vector<vector<int> > Points;
		map<vector<int>,int> PointToIndexMap;
		vector<Edge> Edges;
		vector<Facet> Facets;
		C_Polyhedron CPolyhedron;
		Generator_System LinealitySpace;
		int AffineDimension;
		int SpaceDimension;
};

//------------------------------------------------------------------------------
double GetPolyhedralIntersectionTime();

//------------------------------------------------------------------------------
C_Polyhedron IntersectCones(C_Polyhedron &ph1, C_Polyhedron &ph2);

//------------------------------------------------------------------------------
Cone IntersectCones(Cone C1, Cone C2);

//------------------------------------------------------------------------------
vector<int> GeneratorToPoint(Generator g);

//------------------------------------------------------------------------------
vector<vector<int> > GeneratorSystemToPoints(Generator_System gs);

//------------------------------------------------------------------------------
vector<int> ConstraintToPoint(Constraint c);

//------------------------------------------------------------------------------
Hull NewHull(vector<vector<int> > Points);

//------------------------------------------------------------------------------
vector<Facet> FindFacets(Hull H);

//------------------------------------------------------------------------------
vector<Edge> FindEdges(Hull H);

//------------------------------------------------------------------------------
vector<vector<int> > FindCandidateEdges(Hull H);

//------------------------------------------------------------------------------
int InnerProduct(vector<int> V1, vector<int> V2);

//------------------------------------------------------------------------------
vector<vector<int> > FindInitialForm(vector<vector<int> > &Points, vector<int> &Vector);

//------------------------------------------------------------------------------
C_Polyhedron FindCPolyhedron(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoints(vector<vector<int> > Points);

//------------------------------------------------------------------------------
void PrintPoint(set<int> Point);

//------------------------------------------------------------------------------
void PrintCPolyhedron(C_Polyhedron ph, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
void PrintCPolyhedrons(vector<C_Polyhedron> phs, bool PrintIf0Dim = true);

//------------------------------------------------------------------------------
set<int> IntersectSets(set<int> S1, set<int> S2);

//------------------------------------------------------------------------------
bool SetsDoIntersect(set<int> &S1, set<int> &S2);
