import PolymerLatticeGeneration as poly
from time import *

start = time()

box = poly.Box([100,100,100])

A = 1
B = 2
C = 3

Bnode1 = poly.Node(B)
poly0 = poly.Polymer(Bnode1)
# chain0 = poly.Chain(98,2)
# node0.add_child(chain0)
# chain0.add_child( poly.Node(1) )

for i in range(100):
	if i > 0:
		Bnode1 = poly.Node(B)
		Cnode1.add_child(Bnode1)

	Bnode1.add_child( poly.Node(A) )
	Bnode1.add_child( poly.Node(A) )

	Bnode2 = poly.Node(B)
	Bnode1.add_child(Bnode2)

	Bnode2.add_child(poly.Node(A))
	Bnode2.add_child(poly.Node(A))
	Cnode1 = poly.Node(C)
	Bnode2.add_child(Cnode1)

box.add_polymer(poly0, startPos=[0,0,0])

box.write_box('other-scripts/polymer-Nick_09-25.data')

end = time()
print end-start, 'seconds total'


node0 = poly.Node(1)
poly0 = poly.Polymer(node0)

chain0 = poly.Chain(100, 1)
loop0 = poly.Loop( [chain0] )

node0.add_child( loop0 )



