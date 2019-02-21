import sys
sys.path.insert(0,'../../../../ArbitraryPolymersForLAMMPS')
import PolymerLatticeGeneration as poly
import time


start = time.time()

box = poly.Box([100,100,100])

np=10
node0 = poly.Node(1)
poly0=poly.Polymer(node0)
node0.add_child(poly.Loop([poly.Chain(100,2),poly.Node(1)]))
for _ in range(np):
    box.add_polymer(poly0)

node0 = poly.Node(1)
poly0 = poly.Polymer(node0)
chain0 = poly.Chain(100,2)
node0.add_child(chain0)
chain0.add_child( poly.Node(1) ) 

for _ in range(np):
	box.add_polymer(poly0)


box.write_box('other-scripts/polymer0.data')

print time.time()-start, 'seconds'
