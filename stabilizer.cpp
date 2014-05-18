#include "stabilizer.h"
#include "graphsim.h"
#ifdef WITH_MATHLINK
   #include <mathlink.h>
#endif   
#include <cstdlib>

Stabilizer::Stabilizer (const VertexIndex numQubits_):
   paulis (numQubits_, vector<LocCliffOp> (numQubits_, lco_Id)),
   rowsigns (numQubits_, rp_p1),
   vtxidx (numQubits_, 0)
{
   numQubits = numQubits_;
}


Stabilizer::Stabilizer (const GraphRegister& gr, 
      const hash_set<VertexIndex>& qubits):
   paulis (qubits.size(), vector<LocCliffOp> (qubits.size(), lco_Id)),
   rowsigns (qubits.size()),
   vtxidx (qubits.size())
{
   numQubits = qubits.size ();
   // Build the graph adjacency matrix with Z's and X's in the diagonal
   // and apply the local Clifford unitaries:
   int in = 0;
   for (VtxIdxIter i = qubits.begin(); i != qubits.end(); i++, in++) {
      rowsigns[in] = RightPhase (0);
	   vtxidx[in] = *i;
      int jn = 0;
      for (VtxIdxIter j = qubits.begin(); j != qubits.end(); j++, jn++) {
         if (i==j) {
            paulis[in][jn] = lco_X;
         } else {
           if (gr.vertices[*i].neighbors.find(*j) !=
                 gr.vertices[*i].neighbors.end()) {
              paulis[in][jn] = lco_Z;
           } else {
              paulis[in][jn] = lco_Id;
           }
         }
         // Now the local Clifford unitaries:
         conjugate (in, jn, gr.vertices[*j].byprod);
      }   
   }
}

Stabilizer::Stabilizer (struct QState * qs) :
   paulis (qs->n, vector<LocCliffOp> (qs->n, lco_Id)),
   rowsigns (qs->n),
   vtxidx (qs->n)
{
   static const unsigned char optbl[4] = {0, 1, 3, 2};
   numQubits = qs->n;
   for (int i = 0; i < (int)numQubits; i++) {
      rowsigns[i] = RightPhase (qs->r[numQubits + i]);
	   vtxidx[i] = i;
      for (int j = 0; j < (int)numQubits; j++) {
         bool xhere = ((qs->x [numQubits+i] [j >> 5]) & (1 << (j & 0x1f))) > 0;
         bool zhere = ((qs->z [numQubits+i] [j >> 5]) & (1 << (j & 0x1f))) > 0;
         paulis[i][j] = LocCliffOp (optbl [(zhere<<1) | xhere]);
      }
   }
}


void Stabilizer::add_row (unsigned target, unsigned addend)
{
   //D cerr << "adding row " << addend << " to row " << target << endl;
   for (unsigned col = 0; col < numQubits; col++) {
      rowsigns[target] = rowsigns[target] + 
        LocCliffOp::mult_phase (paulis[target][col], paulis[addend][col]);
      paulis[target][col] = paulis[target][col] * paulis[addend][col];
   }
}

void Stabilizer::conjugate (unsigned row, unsigned col, 
   const LocCliffOp trans)
{
  rowsigns[row] = rowsigns[row] + paulis[row][col].conjugate (trans);   
}

void Stabilizer::conjugate_column (unsigned col, const LocCliffOp trans)
{
   for (unsigned row = 0; row < numQubits; row++) {
      conjugate (row, col, trans);
   }
}

void Stabilizer::print (ostream &os) const
{
   for (unsigned i = 0; i < numQubits; i++) {
      os << rowsigns[i].get_name() << " ";
      for (unsigned j = 0; j < numQubits; j++) {
         os << paulis[i][j].get_name().substr(0,1) << " ";
      }
      os << endl;
   }
}

#ifdef WITH_MATHLINK

MLENV mathlinkenv = 0;
MLINK mathlink = 0;

void close_mathlink (void)
{
   if (mathlink) {
      MLPutFunction (mathlink, "Exit", 0);   
      MLClose (mathlink);
   }
   if (mathlinkenv) {
      MLDeinitialize (mathlinkenv);
   }
}     

void establish_mathlink (void)
{
   long err;
   
   mathlinkenv = MLInitialize (NULL);
   if (!mathlinkenv) {
      cerr << "failed to initialize MathLink library\n";
      exit (1);
   }
   atexit (close_mathlink);
      
   mathlink = MLOpenString (mathlinkenv, "-linkname \"math -mathlink\"", &err);
   // for Windows or Mac, replace line above with:
   // mathlink = MLOpenString (ep, "-linkmode launch", &err);
   
   if (!mathlink) {
      MLDeinitialize (mathlinkenv);
      cerr << "failed to establish MathLink\n";
      exit (1);
   }
   
   MLPutFunction (mathlink, "EvaluatePacket", 1L);
   MLPutFunction (mathlink, "Get", 1L);
   MLPutString (mathlink, "GaussTrans.m");
   MLEndPacket (mathlink);
   
   int pkt;
   while ((pkt = MLNextPacket (mathlink)) && pkt != RETURNPKT)
      MLNewPacket (mathlink);
   if (!pkt) {
      cerr << "no packet.\n";
   }
   MLNewPacket (mathlink);   
   DBGOUT ("MathLink ready.\n");
}   
   

bool Stabilizer::compare (const Stabilizer &othr)
{
   assert (numQubits == othr.numQubits);
   int *thismatr = (int*) malloc (2*numQubits*numQubits * sizeof (int));
   int *othrmatr = (int*) malloc (2*numQubits*numQubits * sizeof (int));
   for (int r = 0; r < numQubits; r++) {
      for (int c = 0; c < numQubits; c++) {
         thismatr [numQubits * 2*r + c] = 
            (paulis[r][c] == lco_X || paulis[r][c] == lco_Y) ? 1 : 0;
         thismatr [numQubits * (2*r+1) + c] = 
            (paulis[r][c] == lco_Z || paulis[r][c] == lco_Y) ? 1 : 0;
         othrmatr [numQubits * 2*r + c] = 
            (othr.paulis[r][c] == lco_X || othr.paulis[r][c] == lco_Y) ? 1 : 0;
         othrmatr [numQubits * (2*r+1) + c] = 
            (othr.paulis[r][c] == lco_Z || othr.paulis[r][c] == lco_Y) ? 1 : 0;
      }
   }
   
   if (!mathlink) {
      establish_mathlink ();
   }
   DBGOUT ("Querying Mathematica\n");
   MLPutFunction (mathlink, "EvaluatePacket", 1L);
   //MLPutFunction (mathlink, "ToString", 1L);
   MLPutFunction (mathlink, "getGaussTrans", 2L);
   long dims[] = {numQubits, 2*numQubits};
   MLPutIntegerArray (mathlink, thismatr, dims, NULL, 2);
   MLPutIntegerArray (mathlink, othrmatr, dims, NULL, 2);
   MLEndPacket (mathlink);
   
   while (MLNextPacket (mathlink) != RETURNPKT) {
      MLNewPacket (mathlink);
      if (MLError (mathlink)) {
         cerr << "MathLink error: " << MLErrorMessage (mathlink) << endl;
         cerr << "Exiting.";
         exit (1);
      }
   }
            
   //char *reply;
   //MLGetString (mathlink, (const char **) &reply);
   //cout << "Mathematika says: " << reply << endl;
   //MLDisownString (mathlink, reply);
   //return true;
   
   int *trans;
   long *rdims;
   char **heads;
   long depth;
   MLGetIntegerArray (mathlink, &trans, &rdims, &heads, &depth);
   DBGOUT ("Got reply.\n");
   assert (depth == 2);
   assert (heads[0] == string ("List"));
   assert (heads[1] == string ("List"));
   assert (rdims[0] == numQubits);
   assert (rdims[1] == numQubits);
   if (trans[0] == -1) {
      DBGOUT ("Stabilizers differ. No transformation found!");
      return false;
   }   
   // What we now have in trans is a description on how to transform
   // this stabilizer to get the other one. trans is to be read as a matrix
   // to be left-multiplied with the stabilizer.
   
   // Do the transformation, i.e. multiply trans with this stab
   // and store the result in trstab.  
   Stabilizer trstab (numQubits);
   for (int r = 0; r < numQubits; r++) {
      for (int c = 0; c < numQubits; c++) {
         assert (trans [r*numQubits + c] == 0 || trans [r*numQubits + c] == 1);
         if (trans [r*numQubits + c]) {
            // add row c of this stab to row r of trstab
            trstab.rowsigns[r] = trstab.rowsigns[r] + rowsigns[c];
            for (int cc = 0; cc < numQubits; cc++) {
               trstab.paulis[r][cc] = trstab.paulis[r][cc] * paulis[c][cc];
               trstab.rowsigns[r] = trstab.rowsigns[r] +
                  LocCliffOp::mult_phase (trstab.paulis[r][cc], paulis[c][cc]);
            }
         }
      }
   }

   // We don't need these anymore:
   free (thismatr);
   free (othrmatr);
   MLDisownIntegerArray (mathlink, trans, rdims, heads, depth);   
   
   //Now, compare the other stab with trstab:
   for (int r = 0; r < numQubits; r++) {
      if (othr.rowsigns[r] != trstab.rowsigns[r]) {
         DBGOUT ("Stabilizer differ: Mismatch of rowsigns in row " << r << endl);
         return false;
      }
      for (int c = 0; c < numQubits; c++) {
         // If this fails there is something wrong with the Mathematica part:
         assert (othr.paulis[r][c] == trstab.paulis[r][c]);
      }
   }
   return true;
}

#else //WITH_MATHLINK 

bool Stabilizer::compare (const Stabilizer &othr)
{
   cerr << "not compiled for use with MathLink";
   exit (1);
}   

#endif //WITH_MATHLINK
