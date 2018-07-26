from matplotlib import pyplot as plt

x1=[]
x2=[]
y1=[]
y2=[]
with open('pairwrite_LJ126.txt', 'r') as f:
	for row in f:
		try:
			x1.append(float(row.split()[1]))
			y1.append(float(row.split()[2]))
		except Exception as e:
			pass
with open('pairwrite_tanh.txt', 'r') as f:
	for row in f:
		try:
			x2.append(float(row.split()[1]))
			y2.append(float(row.split()[2]))
		except Exception as e:
			pass

#print x

plt.plot(x1,y1,'-r', x2,y2,'-b')
plt.title('Pair Write Output from LAMMPS')
plt.ylim(0,1400)
plt.show()