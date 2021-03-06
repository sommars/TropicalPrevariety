#include <ppl.hh>
#include <iostream>
#include <list>
#include <string>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}


class Node {
	public:
		vector<Node> Children;
		bool IsLeaf;
		string Generator;
		Cone LeafCone;
};

//------------------------------------------------------------------------------
void InsertCone(Node &Tree, Cone &LeafCone, int HullIndex);

//------------------------------------------------------------------------------
vector<Cone> GetCones(Node &Tree, int HullIndex);
