import matplotlib.pyplot as plt
import csv
import sys

x = []
y = []

with open('result.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[1]))
        y.append(float(row[0]))

plt.plot(x,y, label = '2D Trajectory' )
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trajectory')
plt.legend()
plt.grid()
plt.show()
