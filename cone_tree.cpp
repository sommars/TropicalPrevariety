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
	for (int i = 0; i != GeneratorStrings.size(); i++) {
		bool FoundCorrectChild = false;
		bool IsLastLoop = i+1 == GeneratorStrings.size();
		vector<Node> *Children = &(*CurrentRoot).Children;
		for (int j = 0; j != (*Children).size(); j++) {
			Node *Child = &((*Children)[j]);
			if ((*Child).Generator == GeneratorStrings[i]) {
				FoundCorrectChild = true;
				CurrentRoot = Child;
				if (IsLastLoop) {
					if ((*Child).IsLeaf) {
						for (size_t k = HullIndex; k != LeafCone.IntersectionIndices.size(); k++) {
							((*Child).LeafCone).IntersectionIndices[k] = IntersectSets(((*Child).LeafCone).IntersectionIndices[k], LeafCone.IntersectionIndices[k]);
						};
					} else {
						(*Child).IsLeaf = true;
						(*Child).LeafCone = LeafCone;
					};
					return;
				};
				break;
			};
		};
		if (!FoundCorrectChild) {
			Node NewNode;
			NewNode.Generator = GeneratorStrings[i];
			// Set the leaf properties correctly, if it is a leaf.
			if (!IsLastLoop) {
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
vector<Cone> GetCones(Node &Tree, int HullIndex) {
	vector<Cone> Result;
	
	if (Tree.IsLeaf) {
		bool ShouldAddCone = true;
		for (size_t k = HullIndex; k != Tree.LeafCone.IntersectionIndices.size(); k++) {
			if (Tree.LeafCone.IntersectionIndices[k].size() == 0) {
				ShouldAddCone = false;
			};
		};
		if (ShouldAddCone) {
			Result.push_back(Tree.LeafCone);
		};
	};
	for (size_t i = 0; i != Tree.Children.size(); i++) {
		vector<Cone> NewCones = GetCones(Tree.Children[i], HullIndex);
		Result.insert(Result.end(), NewCones.begin(), NewCones.end());
	};
	
	return Result;
};
