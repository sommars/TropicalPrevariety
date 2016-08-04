#include <ppl.hh>
#include <iostream>
#include <list>
#include <string>
#include "prevariety_util.h"
#include "polyhedron_tree.h"
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//No node can also be a leaf. If that seems to be the case, then the node that is a leaf but also has a children
//can cut off the children, because it contains its children.
//The number of constraints varies, but it is >=dimension of the initial problem.

void InsertCone(Node &Tree, Cone LeafCone) {
	vector<string> ConstraintStrings;
	// Convert the constraints to a list of strings. Sort them by lex.
};

//------------------------------------------------------------------------------
vector<Cone> GetCones(Node &Tree) {
	vector<Cone> Result;
	
	if (Tree.IsLeaf) {
		Result.push_back(Tree.LeafCone);
	} else {
		for (int i = 0; i != Tree.Children.size(); i++) {
			vector<Cone> NewCones = GetCones(Tree.Children[i]);
			Result.insert(Result.end(), NewCones.begin(), NewCones.end());
		};
	};
	
	return Result;
};
