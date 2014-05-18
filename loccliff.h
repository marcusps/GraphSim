/*!\file
// loccliff.h -- a class for operators in the Clifford group

// version v0.11, of 2005-01-27

// Copyright (C) 2004  Simon Anders  <sanders@fs.tum.de>
// Institute of Theoretical Physics, University of Innsbruck, Austria

// ----------

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

// ----------
*/

#ifndef LOCCLIFF_H
#define LOCCLIFF_H

#include <string>
#include <vector>

#include <cassert>

using namespace std;


// ----=== RightPhase ===----

//! A phase, which is 0, 90, 180, or 270 degrees.
struct RightPhase {
   
   /*! to be read modulo 4, i.e. & 0x03
   The field 'ph' indicated the angle with a value 0, 1, 2, or 3,
   corresponding to 0, 90, 180, or 270 degrees
   Do not rely on ph<3, always read it out anded with 0x03. */
   unsigned short ph;

   RightPhase (void);               //!< initializes as 0 degrees
   RightPhase (unsigned short ph_); //!< takes 0..3 as argument and writes it to ph
   
   string get_name (void) const;    //!< returns " +", " i", " -", or "-i"
};


// ----=== RightMatrix ===----

#ifndef SWIG
struct RightMatrix {
   bool sqrt2norm;
   bool ampls[2][2];
   RightPhase phases[2][2];
   bool apply_on_state (vector<bool>::reference ampl1, 
      vector<bool>::reference ampl2, RightPhase& ph1, RightPhase& ph2);
};
#endif


// ----=== LocCliffOp ===----

//! An operator in the local Clifford group.
struct LocCliffOp {

  //! The field 'op' identifies the operator. 0 is identity I, 1 is Pauli X,
  //! 2 is Pauli Y, 3 is Pauli Z, 4 is I^B, 5 is I^C etc.
  unsigned short op;

  //! constructor, takes an integer in 0..23.
  LocCliffOp (unsigned short int op_);

  //! constructor, takes a sign symbol in 0..3 (for I, X, Y, Z) and a
  /*! permutation symbol 0..5 (for A, B, ..., F) */
  LocCliffOp (unsigned short int signsymb, unsigned short int permsymb);
  
  //! returns something like e.g. "YC" for Hadamard=Y^C
  string get_name (void) const;
  
  //! replaces op by trans * op * trans^dagger and returns a phase,
  /*! either +1 or -1 (as RightPhase(0) or RightPhase(2)) */
  RightPhase conjugate (const LocCliffOp trans);
  
  //! returns the Hermitian adjoint of op
  LocCliffOp herm_adjoint (void) const;
  
  //! returns the phase of the multiplication of op1 * op2
  static RightPhase mult_phase (LocCliffOp op1, LocCliffOp op2);
  
  //! returns True if op is XA or YA.
  bool isXY (void) const;
  
  //! returns True if op is an operator diagonal in the computational basis
  bool is_diagonal (void) const;
  
  RightMatrix get_matrix (void) const;
};


// ----=== operators ===----

#ifndef SWIG   // I should add names for the operators so that SWIG
               // can handle them

//! muliplies two local Clifford operators. Not commutative!
LocCliffOp operator* (LocCliffOp a, LocCliffOp b);

//! Check two LC operators for equality
bool operator== (LocCliffOp a, LocCliffOp b);
//! Check two LC operators for inequality
bool operator!= (LocCliffOp a, LocCliffOp b);

//! Adds two phases (modulo 4).
RightPhase operator+ (RightPhase ph1, RightPhase ph2);

//! Check two phases for equality
bool operator== (RightPhase ph1, RightPhase ph2);
//! Check two phases for inequality
bool operator!= (RightPhase ph1, RightPhase ph2);

#endif //SWIG

// ---=== constants ===---

//! There are 24 local Clifford operators.
const unsigned short num_LocCliffOps = 24;

// special local Clifford operators:
const LocCliffOp lco_Id   = 0;        //!< Identity
const LocCliffOp lco_X    = 1;        //!< Pauli X
const LocCliffOp lco_Y    = 2;        //!< Pauli Y 
const LocCliffOp lco_Z    = 3;        //!< Pauli Z
const LocCliffOp lco_H    = 10;       //!< Hadamard
const LocCliffOp lco_spiZ = 5;        //!< Sqrt (+iZ)
const LocCliffOp lco_smiZ = 6;        //!< Sqrt (-iZ)
const LocCliffOp lco_spiY = 11;       //!< Sqrt (+iY)
const LocCliffOp lco_smiY = 9;        //!< Sqrt (-iY)
const LocCliffOp lco_spiX = 14;       //!< Sqrt (+iX)
const LocCliffOp lco_smiX = 15;       //!< Sqrt (-iX)
const LocCliffOp lco_S    = lco_smiZ; //!< Pi/4 phase rot
const LocCliffOp lco_Sh   = lco_spiZ; //!< herm. conj. of Pi/4 phase rot

// the right angles:
const RightPhase rp_p1 = 0;         //!< +1
const RightPhase rp_pI = 1;         //!< +i
const RightPhase rp_m1 = 2;         //!< -1
const RightPhase rp_mI = 3;         //!< -i


// ---=== internals ===---


#ifndef SWIG

//! Inline functions and tables needed by them. Consider them as private.
// (all others functions in loccliff.cpp)
namespace loccliff_tables {

const unsigned short meas_conj_tbl [3] [num_LocCliffOps] = 
{{1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3}, 
 {2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1}, 
 {3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2}};

const unsigned short lco_mult_tbl [num_LocCliffOps] [num_LocCliffOps] = 
#include "multtbl.tbl"
;    

const short adj_tbl [24] =
   { 0,  1,  2,  3,
     4,  6,  5,  7,
     8, 11, 10,  9,
    12, 13, 15, 14,
    20, 22, 23, 21,
    16, 19, 17, 18};    
    
 const short phase_tbl[4][4] =
  {{0, 0, 0, 0},
   {0, 0, 1, 3},
   {0, 3, 0, 1},
   {0, 1, 3, 0}};

} // end namespace loccliff_tables


inline LocCliffOp::LocCliffOp (unsigned short int op_) 
{
   op = op_;
}

inline LocCliffOp::LocCliffOp (unsigned short int signsymb, 
      unsigned short int permsymb) 
{
   assert (signsymb < 4 && permsymb < 6);
   op = permsymb * 4 + signsymb;
}
 
inline LocCliffOp LocCliffOp::herm_adjoint (void) const
{
   return loccliff_tables::adj_tbl [op];
}

inline RightPhase LocCliffOp::mult_phase (LocCliffOp op1, LocCliffOp op2)
{
   assert (op1.op <= lco_Z.op && op2.op <= lco_Z.op);
   return RightPhase (loccliff_tables::phase_tbl[op1.op][op2.op]);
}

//!
inline bool LocCliffOp::isXY (void) const
{
   return (op == lco_X.op || op == lco_Y.op);
}

inline bool LocCliffOp::is_diagonal (void) const
{
   return (op == lco_Id.op   || op == lco_Z.op || 
           op == lco_smiZ.op || op == lco_spiZ.op);
}

inline LocCliffOp operator* (LocCliffOp a, LocCliffOp b) 
{ 
  return LocCliffOp (loccliff_tables::lco_mult_tbl [a.op] [b.op]);
}

inline bool operator== (LocCliffOp a, LocCliffOp b) 
{
  return a.op == b.op;
}

inline bool operator!= (LocCliffOp a, LocCliffOp b) 
{
  return a.op != b.op;
}

inline RightPhase::RightPhase (void) {
  ph = 0;
}

inline RightPhase::RightPhase (unsigned short ph_) {
  ph = ph_;
}

inline RightPhase operator+ (RightPhase ph1, RightPhase ph2)
{
   return RightPhase ((ph1.ph + ph2.ph) & 0x03);
}

inline bool operator== (RightPhase a, RightPhase b) 
{
  return ((a.ph ^ b.ph) & 0x03) == 0;
}
inline bool operator!= (RightPhase a, RightPhase b) 
{
  return ((a.ph ^ b.ph) & 0x03) != 0;
}
#endif //SWIG

#endif //LOCCLIFF_H
