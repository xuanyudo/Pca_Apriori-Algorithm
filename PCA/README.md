# There are two version of implementations.
## pca1.py
pca1.py is capable to run PCA dimension algorithm on dataset you pass in with argument. Usage: <br/>
#### python3 pca1.py [filename] <br/>
Where as [filename] is your path to your input file with the same format as pca_a.txt, pca_b.txt or pca_c.txt.

## pca_with_numpy.py
pca1.py is capable to run PCA, SVD or TSNE algorithm on dataset you pass in with argument. Usage: <br/>
#### python3 pca_with_numpy.py [filename] [pca|svd|tsne]<br/>
Where as [filename] is your path to your input file with the same format as pca_a.txt, pca_b.txt or pca_c.txt.
and [pca|svd|tsne] is three algorithms you can pick from. 
Omitting this option will result in running pca algorithm in default, and inputting wrong option will result an error.

## Result
Program will result in generating a graph after reduction process. and the graph will be show using matplotlib library.
