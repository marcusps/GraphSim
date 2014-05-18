// graphsim.h

/*!\mainpage
This is 'graphsim', a simulator for stabilizer quantum cuircuits using the
graph state formalism.

graphsim version 0.10, dated 2005-Jan-27

(c) Simon Anders (sanders@fs.tum.de), University of Innsbruck

See the article

  S. Anders, H. J. Briegel: 
  Fast simulation of Stabilizer Circuits using a Graph States Formalism
  quant-ph/0504117

for a description of this software.

----

// This file is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.

// This file is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this file; see the file COPYING.  If not, browse to 
// http://www.fsf.org/licenses/gpl.html

*/

/*!\file
This header file defines the main interface of graphsim.

(c) Simon Anders, University of Innsbruck, 2005
released under GPL.

Note: If you have trouble compiling this, please note:
This file uses the "hash_set" template, which is an extension
to the Standard C++ Library and the Standard Template Library,
specified by SGI. Most C++ compilers have this template.
For GNU C++, the header file <ext/hash_set> hasa to be included
and hash_set has to be prefixed with namespace __gnu_cxx.
If you use another compiler, you might have to change the include
file and the namespace identifier.
*/

#ifndef GRAPHSIM_H
#define GRAPHSIM_H

//The following directives are for SWIG, see http://www.swig.org
#ifdef SWIG
%module (docstring="Graph State Stabilizer Simulator -- S. Anders") graphsim
%{
#include "graphsim.h"
#include "loccliff.h"
#include "stabilizer.h"
%}
%include "cpointer.i"
%pointer_class (bool, boolpc);
%feature ("autodoc", "1");
%rename (print_tbl) Stabilizer::print;
%rename (print_kets) CBDecomposition::print;
%include "graphsim.h"
%include "loccliff.h"
%include "stabilizer.h"
#endif // SWIG

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>

#include "loccliff.h"
#include "stabilizer.h"

#include <cstring>
#include <cstdlib>

using namespace std;

/*! All vertices in a graph state are numbered beginning with 0. To specify
auch an index, the type VertexIndex (which is just unsigned long) is always
used */
typedef unsigned long VertexIndex;

/*! A GraphRegister object maintains a list of its vertices (qubits), each 
described by an object of the class QubitVertex described here.*/
struct QubitVertex {
   /*!byprod is the vertex operator (VOp) associated with the qubit (the name 
   stems from the term 'byproduct operator' used for the similar concept in 
   the one-way quantum computer.*/
   LocCliffOp byprod;
   /*! neigbors is the adjacency list for this vertex */
   hash_set<VertexIndex> neighbors;
   /*! Upon construction, a qubit vertex is initialised with the Hadamard
   operation as VOp, and with wmpty neighbor list. This makes it represent
   a |0>. */
   QubitVertex (void)
     : byprod (lco_H) {};
};
   
#ifndef SWIG
/*! This structure is only for internal use for the cphase functions. */
struct ConnectionInfo {
   bool wasEdge;
   bool non1;
   bool non2;
};
#endif

//! A quantum register.
/*! GraphRegister is the central class of graphsim. It represents a register of qubits
that can be entangled with each other. It offers functions to initialize the register,
let gates operate on the qubits, do measurements and print out the state. */
class GraphRegister {
  public:
   /*! This vector stores all the qubits, represented as QubitVertex objects. The index
   of the vector is usually taken as of type VertexIndex. */
   vector<QubitVertex> vertices;
   GraphRegister (VertexIndex numQubits, int randomize = -1);
   GraphRegister (GraphRegister& gr);
   ~GraphRegister () {};
   void local_op (VertexIndex v, LocCliffOp o);
   void hadamard (VertexIndex v);
   void phaserot (VertexIndex v);
   void bitflip (VertexIndex v);
   void phaseflip (VertexIndex v);
   void cphase (VertexIndex v1, VertexIndex v2);      
   void cnot (VertexIndex vc, VertexIndex vt);
   int measure (VertexIndex v, LocCliffOp basis = lco_Z, 
      bool* determined = NULL, int force = -1);
   Stabilizer & get_full_stabilizer (void) const;
   void invert_neighborhood (VertexIndex v);
   void print_adj_list (ostream& os = cout) const;
   void print_adj_list_line (ostream& os, VertexIndex i) const;
   void print_stabilizer (ostream& os = cout) const;
  private:
   void add_edge (VertexIndex v1, VertexIndex v2);
   void del_edge (VertexIndex v1, VertexIndex v2);
   void toggle_edge (VertexIndex v1, VertexIndex v2);
   int graph_Z_measure (VertexIndex v, int force = -1);
   int graph_Y_measure (VertexIndex v, int force = -1);
   int graph_X_measure (VertexIndex v, bool* determined = NULL, int force = -1);
   void toggle_edges (const hash_set<VertexIndex> vs1, 
      const hash_set<VertexIndex> vs2);      
   bool remove_byprod_op (VertexIndex v, VertexIndex use_not);
   void cphase_with_table (VertexIndex v1, VertexIndex v2);
   ConnectionInfo getConnectionInfo (VertexIndex v1, VertexIndex v2);  
};


#ifndef SWIG
/*! As we often iterate over sublists of GraphRegister::vertices, this
iterator typedef is a handy abbreviation. */
typedef vector<QubitVertex>::iterator VertexIter;
/*! Another iterator, this one for the adjacency lists QubitVertex::neigbors,
and subsets. */
typedef hash_set<VertexIndex>::iterator VtxIdxIter;

/*! A constant version of VertexIter */
typedef vector<QubitVertex>::const_iterator VertexIterConst;
/*! A constant version of VtxIdxIter */
typedef hash_set<VertexIndex>::const_iterator VtxIdxIterConst;


/*! Apply the local (i.e. single-qubit) operation o on vertex v. */
inline void GraphRegister::local_op (VertexIndex v, LocCliffOp o) {
   vertices[v].byprod = o * vertices[v].byprod;
}

/*! Apply a Hadamard gate on vertex v */
inline void GraphRegister::hadamard (VertexIndex v) {
   local_op (v, lco_H);
} 

/*! Apply a phaserot gate on vertex v. Phaserot means the gate
S = |0><0| + i |1><1|. */
inline void GraphRegister::phaserot (VertexIndex v) {
   local_op (v, lco_S);
} 

/*! Apply a bitflip gate (i.e. a Pauli X) on vertex v */
inline void GraphRegister::bitflip (VertexIndex v) {
   local_op (v, lco_X);
} 

/*! Apply a phaseflip gate (i.e. a Pauli Z) on vertex v */
inline void GraphRegister::phaseflip (VertexIndex v) {
   local_op (v, lco_Z);
} 

#endif //SWIG

//#define DEBUGOUTPUT

#ifdef DEBUGOUTPUT
   #define DBGOUT(a) cout << a
#else
   #define DBGOUT(a) 
#endif



#endif //GRAPHSIM_H
