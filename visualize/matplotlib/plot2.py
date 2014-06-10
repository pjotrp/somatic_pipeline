import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import csv

# image = np.random.uniform(size=(10, 10))
data = [[1,0,2],[0,2,2],[0,3,4]]
# ax.set_title('Mutated genes')

column_labels = ['g1','g2','g3']
row_labels = ['p1','p2','p3']

file_name = 'data/test.csv'
# data = np.loadtxt(open(file_name), delimiter=',')
reader = csv.reader(open('data/test.csv'),delimiter=',')
x=list(reader)
matrix = np.array(x)
rownames = [row[0] for row in matrix[1:]]
colnames = matrix[0][1:]

print matrix
data = np.delete(matrix,0,0)
print data
data = np.delete(data,0,1)
print data
data[data=='']='0'

# print row_labels, column_labels
data=data.astype(int)
print data
# print data

fig, ax = plt.subplots()

ax.imshow(data, cmap=matplotlib.cm.Blues, interpolation='nearest')
# Move left and bottom spines outward by 10 points
# ax.spines['left'].set_position(('outward', 3))
# ax.spines['bottom'].set_position(('outward', 3))
# Hide the right and top spines
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
# ax.yaxis.set_ticks_position('left')
# ax.xaxis.set_ticks_position('top')

# want a more natural, table-like display
# ax.invert_yaxis()
ax.xaxis.tick_top()

print rownames
ax.set_xticks(np.arange(len(colnames)), minor=False)
ax.set_yticks(np.arange(len(rownames)), minor=False)
ax.set_xticklabels(colnames, minor=False)
ax.set_yticklabels(rownames, minor=False)

# ax.set_title("Mains power stability")    
ax.set_xlabel('Samples')
ax.set_ylabel('Genes')

# ax.autoscale(tight=True)
# ax.grid()

# plt.show()
# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(14.5,2.5)
plt.savefig('test2png.png',dpi=100)

plt.pause(0.1) # <-------
raw_input("<Hit Enter To Close>")
plt.close(fig)
