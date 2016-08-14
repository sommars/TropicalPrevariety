#include <ppl.hh>
#include <iostream>
#include <list>
#include <string>
#include "prevariety_util.h"
#include "cone_tree.h"
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

//No node can also be a leaf. If that seems to be the case, then the node that is a leaf but also has a children
//can cut off the children, because it contains its children.
//The number of constraints varies, but it is >=dimension of the initial problem.

//------------------------------------------------------------------------------
void InsertCone(Node &Tree, Cone &LeafCone, int HullIndex) {
	vector<string> GeneratorStrings;
	// Convert the generators to a list of strings. Sort them by lex
	Generator_System cs = LeafCone.Polyhedron.minimized_generators();
	for (Generator_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		stringstream ss;
		ss << (*i);
		GeneratorStrings.push_back(ss.str());
	};
	sort(GeneratorStrings.begin(),GeneratorStrings.end());
	
	// Bite off p(0).
	GeneratorStrings.erase(remove(GeneratorStrings.begin(), GeneratorStrings.end(), "p(0)"),GeneratorStrings.end());
	Node *CurrentRoot = &Tree;
	for (size_t i = 0; i != GeneratorStrings.size(); i++) {
		Node Child;
		Child.Generator = GeneratorStrings[i];
		if (i+1 == GeneratorStrings.size()) {
			Child.IsLeaf = true;
			Child.LeafCone = LeafCone;
		} else {
			Child.IsLeaf = false;
		};
		pair<map<string, Node>::iterator,bool> ret;
		ret = (*CurrentRoot).Children.insert(pair<string,Node>(GeneratorStrings[i], Child));
		map<string, Node> *Children = &(*CurrentRoot).Children;
		
		CurrentRoot = ((*Children)[GeneratorStrings[i]]);
	};
};

//------------------------------------------------------------------------------
vector<Cone> GetCones(Node &Tree) {
	//TODO: CHANGE GETCONES: make the list at the same time as the tree is made. can't store index because that will get stale. do erase find paradigm.

	vector<Cone> Result;
	/*
	if (Tree.IsLeaf) {
		Result.push_back(Tree.LeafCone);
	};
	set<Node>::iterator ChildItr;
	for (ChildItr = Tree.Children.begin(); ChildItr != Tree.Children.end(); ChildItr++) {
		Node Child = *ChildItr;
		vector<Cone> NewCones = GetCones(Child);
		Result.insert(Result.end(), NewCones.begin(), NewCones.end());
	};*/
	return Result;
};
