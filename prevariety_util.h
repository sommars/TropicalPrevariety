#include <ppl.hh>
#include <iostream>
#include <list>
using namespace std;
using namespace Parma_Polyhedra_Library;
namespace Parma_Polyhedra_Library {using IO_Operators::operator<<;}

int InnerProduct(list<int> V1, list<int> V2);

list<list<int> > FindInitialForm(list<list<int> > Points, list<int> Vector);

C_Polyhedron FindCPolyhedron(list<list<int> > Points);
