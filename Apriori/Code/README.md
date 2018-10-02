## apriori.py
#### Highly suggest using this one instead of another version since it is much faster and flexible. Usage:
#### python3 apriori.py [filename] [support] [confidence]

Note that you need all three options in order to run <br/>
[filename] indicates the dataset you want to pass in to the algorithm. please format the dataset same way as associationruletestdata.txt.<br/>
[support] numeric value of how much support. support should be set in between [0,1]<br/>
[confidence] numeric value of how much support. confidence should be set in between [0,1]<br/>

Go to line 119 and change the template accordingly. use self.template[num] and template rules to run template. Note that for search items it should be all upper case. 
Or other way, you can use ap.template[num] below ap.apriori(sup, conf), it would return the same (result, count) where as result is a list of rules satisfy the template and count is total number of them. Note that for search items it should be all upper case. 
