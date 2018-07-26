from matplotlib import pyplot as plt

x =[]
y=[]
with open('pair11.txt', 'r') as f:
	for row in f:
		try:
			x.append(float(row.split()[1]))
			y.append(float(row.split()[2]))
		except Exception as e:
			pass


#print x

plt.plot(x,y)
plt.show()