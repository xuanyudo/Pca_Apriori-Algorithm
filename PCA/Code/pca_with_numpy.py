import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import sys
from numpy import genfromtxt

from sklearn.manifold import TSNE


dimensions = 2
testing_sets = np.array([[19, 63], [39, 74], [30, 87], [30, 23],[15, 35], [15, 43], [15, 32], [30, 73]])
	
def PCA(data):
	# calculate set with mean corrected
	corrected_sets = data - data.mean(0)

	# calcialte covariance matrix
	cov_matrix = np.cov(corrected_sets.transpose(), bias=True)
	# get eigens
	eigval, eigvec = np.linalg.eig(cov_matrix)

	# sory by eigen value
	idx = eigval.argsort()[::-1]   
	eigval = eigval[idx]
	eigvec = eigvec[:,idx]
	print(eigvec)
	# multiply to our matrix
	y0 = np.sum(corrected_sets*eigvec[0], axis=1)
	y1 = np.sum(corrected_sets*eigvec[1], axis=1)


	return y0, y1
	
def plot(y0, y1, labels):
	df = pd.DataFrame(dict(x=y0, y=y1, label=labels))
	groups = df.groupby('label')
	fig, ax = plt.subplots()
	ax.margins(0.05)
	for name, group in groups:
	    ax.plot(group.x, group.y, marker='o', linestyle='', ms=4, label=name.decode("utf-8"))
	ax.legend(numpoints=1, loc='upper right')

	plt.show()


def tsneformat(data):
	y = TSNE(n_components=2).fit_transform(formattedData)
	return y[:,0], y[:, 1]

filename = sys.argv[1]

if len(sys.argv)<3:
	algorithm = 'pca'
else:
	algorithm = sys.argv[2]
my_data = genfromtxt(filename, delimiter='\t', dtype="|S10")

formattedData = my_data[:, :-1:].astype(np.float)

if algorithm == 'pca':
	y0, y1 = PCA(formattedData)

elif algorithm == 'svd':
	u, s, vh = np.linalg.svd(my_data[:, :-1:].astype(np.float))
	y0, y1 = u[:,0], u[:,1]
elif algorithm == 'tsne':
	y0, y1 = tsneformat(formattedData)

labels = my_data[:, -1::].transpose()[0]

plot(y0, y1, labels)

