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

void InsertCone(Node &Tree, Cone &LeafCone, int HullIndex) {
	vector<string> ConstraintStrings;
	// Convert the constraints to a list of strings. Sort them by lex
	Constraint_System cs = LeafCone.Polyhedron.minimized_constraints();
	for (Constraint_System::const_iterator i = cs.begin(),
	cs_end = cs.end(); i != cs_end; ++i) {
		stringstream ss;
		ss << (*i);
		ConstraintStrings.push_back(ss.str());
	};
	sort(ConstraintStrings.begin(),ConstraintStrings.end());
	Node *CurrentRoot = &Tree;
	for (int i = 0; i != ConstraintStrings.size(); i++) {
		bool FoundCorrectChild = false;
		vector<Node> *Children = &(*CurrentRoot).Children;
		for (int j = 0; j != (*Children).size(); j++) {
			Node *Child = &((*Children)[j]);
			if ((*Child).Constraint == ConstraintStrings[i]) {
				FoundCorrectChild = true;
				CurrentRoot = &(*Child);
				if ((*Child).IsLeaf) {
					// Either the cone we are trying to add is the same cone, or it's smaller.
					// If it's smaller, we don't want to add it.
					if (i+1 == ConstraintStrings.size()) {
						// Intersect all of the pre-intersect indices.
						for (size_t k = HullIndex; k != LeafCone.IntersectionIndices.size(); k++) {
							((*Child).LeafCone).IntersectionIndices[k] = IntersectSets(((*Child).LeafCone).IntersectionIndices[k], LeafCone.IntersectionIndices[k]);
						};
					};
					return;
				} else {
					if (i+1 == ConstraintStrings.size()) {
						(*Child).IsLeaf = true;
						(*Child).LeafCone = LeafCone;
						(*Child).Children.clear();
						return;
					};
				};
				break;
			};
		};
		if (!FoundCorrectChild) {
			Node NewNode;
			NewNode.Constraint = ConstraintStrings[i];
			// Set the leaf properties correctly, if it is a leaf.
			if (i+1 != ConstraintStrings.size()) {
				NewNode.IsLeaf = false;
			} else {
				NewNode.IsLeaf = true;
				NewNode.LeafCone = LeafCone;
			};
			(*CurrentRoot).Children.push_back(NewNode);
			CurrentRoot = &((*CurrentRoot).Children[(*CurrentRoot).Children.size()-1]);
		};
	};
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
