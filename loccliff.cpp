// loccliff.cpp

// version v0.10, of 2004-10-28

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


#include "loccliff.h"
#include <iostream>

using namespace std;

string LocCliffOp::get_name (void) const
{
   static const char* paulinames[] = {"I", "X", "Y", "Z"};
   return string (paulinames[op & 0x03]) + (char) ('A' + op / 4);
}

RightPhase LocCliffOp::conjugate (const LocCliffOp trans) {
  //If *this is the identity, we don't have to do anything
  if (*this == lco_Id) {
     return RightPhase (0);
  }
  //This is meant to be used only if *this is a Pauli:
  assert (op >= lco_X.op && op <= lco_Z.op);
  // First the sign:
  RightPhase zeta;
  if ((trans.op & 0x03) == 0 || (trans.op & 0x03) == op) {
     // zeta = + sgn pi
     // sgn pi = -1 iff trans.op >= 4 && trans.op <= 15
     if (trans.op >= 4 && trans.op <= 15) {
       zeta = RightPhase (2);
     } else {
       zeta = RightPhase (0);
     }
  } else {
     // zeta = - sgn pi
     // sgn pi = -1 iff trans.op >= 4 && trans.op <= 15
     if (trans.op >= 4 && trans.op <= 15) {
       zeta = RightPhase (0);
     } else {
       zeta = RightPhase (2);
     }
  }
  // Now the operator:
  // First check the table (to be removed!):
  assert (loccliff_tables::meas_conj_tbl [op-lco_X.op] [trans.op] 
    == trans * op * trans.herm_adjoint()); 
  op = loccliff_tables::meas_conj_tbl [op-lco_X.op] [trans.op];
  return zeta;
}

RightMatrix LocCliffOp::get_matrix (void) const
{
   const short matrices[24][2][2] = 
    {{{0, -1}, {-1, 0}}, {{-1, 0}, {0, -1}}, {{-1, 3}, {1, -1}}, {{0, -1}, {-1, 
      2}}, {{-1, 0}, {3, -1}}, {{0, -1}, {-1, 3}}, {{0, -1}, {-1, 1}}, {{-1, 
      0}, {1, -1}}, {{0, 2}, {2, 2}}, {{0, 2}, {0, 0}}, {{0, 0}, {0, 2}}, {{0,
       0}, {2, 0}}, {{0, 1}, {3, 2}}, {{0, 3}, {1, 2}}, {{0, 1}, {1, 0}}, {{0,
       3}, {3, 0}}, {{0, 3}, {0, 1}}, {{0, 1}, {2, 1}}, {{0, 3}, {2, 3}}, {{0,
       1}, {0, 3}}, {{0, 0}, {1, 3}}, {{0, 0}, {3, 1}}, {{0, 2}, {3, 3}}, {{0,
       2}, {1, 1}}};
   RightMatrix rm;
   rm.sqrt2norm = true;
   for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
         if (matrices[op][i][j] == -1) {
            rm.ampls[i][j] = false;
            rm.sqrt2norm = false;
            rm.phases[i][j] = RightPhase (0);
         } else {
            rm.ampls[i][j] = true;
            rm.phases[i][j] = RightPhase (matrices[op][i][j]);
         }
      }
   }
   return rm;
}


string RightPhase::get_name (void) const
{
   const char* names[] = {"  ", " i", " -", "-i"};
   return string (names[ph & 0x03]);
}



bool RightMatrix::apply_on_state (vector<bool>::reference ampl1, 
      vector<bool>::reference ampl2, RightPhase& ph1, RightPhase& ph2)
{
   vector<bool>::reference amplV[2] = {ampl1, ampl2};
   RightPhase *phV[2] = {&ph1, &ph2};
   RightMatrix sum;
   bool diag[2] = {false, false};
   
   for (int r = 0; r < 2; r++) {
      for (int c = 0; c < 2; c++) {
         sum.ampls[r][c] = ampls[r][c] && amplV[c];
         sum.phases[r][c] = phases[r][c] + *phV[c];
      }
   }
   for (int r = 0; r < 2; r++) {
      amplV[r] = sum.ampls[r][0] || sum.ampls[r][1];
      if (!amplV[r]) {
         continue;
      }
      if (! (sum.ampls[r][0] && sum.ampls[r][1])) {
         // not both ampls present -> just copy from one
         *phV[r] = sum.phases[r][sum.ampls[r][0] ? 0: 1];
      } else {
         // both ampls present. We have to add them
         switch ((sum.phases[r][1].ph - sum.phases[r][0].ph + 4) % 4) {
            case 0:
               // They are the same. Take one:
               *phV[r] = sum.phases[r][0]; 
               break;
            case 2:
               // They cancel
               amplV[r] = false;
               break;
            case 1: 
               *phV[r] = sum.phases[r][0]; 
               diag[r] = true;
               break;
            case 3: 
               *phV[r] = sum.phases[r][1]; 
               diag[r] = true;
               break;
            default: assert (0);                  
         }
      }
   }
   if (amplV[0] && amplV[1]) {
      assert (diag[0] == diag[1]);
   }
   assert (amplV[0] || amplV[1]);
   return amplV[0] ? diag[0] : diag[1];  
}
