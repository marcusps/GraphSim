#include <iostream>
#include <string>

#include "graphsim.h"

using namespace std;

const int nbr_of_qubits = 200;

int main (int, char**)
{
   //srandom (time (NULL));
   GraphRegister gr (nbr_of_qubits);
   
   for (int iter = 0; iter < 1500; iter++) {
      VertexIndex qubit = random () / (RAND_MAX / nbr_of_qubits);
      int what = random () / (RAND_MAX / 15);
      switch (what) {
         case 0: gr.hadamard (qubit); break;
         case 1: gr.local_op (qubit, lco_S); break;
         case 2: gr.local_op (qubit, lco_Z); break;
         case 3: gr.local_op (qubit, lco_X); break;
         default:
            VertexIndex qubit2 = random () / (RAND_MAX / nbr_of_qubits);
            if (qubit2 == qubit) {
               continue;
            }
            if (random() > RAND_MAX/2) {
               gr.cphase (qubit, qubit2);
            } else {
               gr.cnot (qubit, qubit2);
            }
      }
      if (iter % 100 == 0) {
         cout << iter << endl;
      }
   }
   cout << endl;
   gr.print_adj_list ();
}
