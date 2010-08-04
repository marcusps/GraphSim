import graphsim
import random

# we need a quantum register with 7 qubits:
gr = graphsim.GraphRegister (8)

gr.hadamard (4)
gr.hadamard (5)
gr.hadamard (6)
gr.cnot (6, 3)
gr.cnot (6, 1)
gr.cnot (6, 0)
gr.cnot (5, 3)
gr.cnot (5, 2)
gr.cnot (5, 0)
gr.cnot (4, 3)
gr.cnot (4, 2)
gr.cnot (4, 1)

for i in xrange (7):
   gr.cnot (i, 7)

print gr.measure (7)

gr.print_adj_list ()
gr.print_stabilizer ()
   
