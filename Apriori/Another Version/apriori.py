
from numpy import genfromtxt
import numpy as np
from functools import reduce
import hashlib
import pickle

support = 0.3
confidence = 0.7
testFileName = 'associationruletestdata'
#testFileName = 'test'
testFile = testFileName+'.txt'
saveFileName = "savedsupport_"+str(support)+"confid_"+str(confidence)+testFileName
freqSetFileName = "saved/"+saveFileName+".pickle"
ruleFileName = "saved/"+saveFileName+"rule.picklerule";

freqPrintFileName = saveFileName+".txt"
rulePrintFileName = saveFileName+"rule.txt";
ruleTemplatePrintFileName = saveFileName+"Templates.txt";

# input: recursive set, data for searching, support percentage
# use recursive breadth first search approach
def aporoi_freq_set(seto, formatted_data, support):
	if(len(seto) == 0):
		return seto
	newSet = set()
	for setobject in seto:
		bary, allItemNames = setobject.return_vars()
		for bstrName in allItemNames:
			aary = set(bary)
			aary.add(bstrName)
			freq, count = get_if_freq(formatted_data, aary, support)
			setRemain = set(allItemNames)
			setRemain.remove(bstrName)
			if(freq):
				newSet.add(dataSet(aary, setRemain, count))
	
	return newSet.union(aporoi_freq_set(newSet, formatted_data, support))

# function reason why this whole program is slow
def get_if_freq(formatted_data, bstrset, support):
	total = len(formatted_data)
	cum = 0
	for colu in formatted_data:
		if bstrset.issubset(colu):
			cum+=1
	return cum/total >= support, cum


def rule_gen(freqTable, confidence):
	rules = set();
	for rule in freqTable:
		#print("#########")
		#print(rule.return_vars())
		#print("#########")
		rules = rules.union(aporoi_single_rule_gen(set([rule]), freqTable, confidence))
	print(len(rules))
	return rules

def aporoi_single_rule_gen(ruleSet, ruleTable, confidence):
	if(len(ruleSet) == 0):
		return ruleSet
	newSet = set()
	for ruleobject in ruleSet:
		body, head = ruleobject.return_vars()
		if(len(head) == 1):
			continue
		therule = ruleobject.return_rule()
		for bstrName in head:
			aset = set(body)
			aset.add(bstrName)
			headset = set(head)
			headset.remove(bstrName)
			conf = get_if_conf(ruleTable, therule, headset, confidence)
			if(conf):
				newSet.add(rule(headset, aset))
		
	return newSet.union(aporoi_single_rule_gen(newSet, ruleTable, confidence))

def get_if_conf(ruleTable, therule, headset, confidence):
	#print(therule)
	countRule = ruleTable[rule(therule, set())].return_count();
	countHead = ruleTable[rule(headset, set())].return_count();
	assert countRule <= countHead
	return countRule/countHead >= confidence

def formatData(dataset):
	for index, bstr in np.ndenumerate(dataset):
		if(bstr == b"Up" or bstr == b"Down"):
			dataset[index] = b"G%d_%s" % (index[1]+1, bstr)
	return dataset

# turn freq set into a table
def formatFreqSet(freqset):
	table = {}
	for data in freqset:
		array, rest = data.return_vars()
		table[rule(set(array), set())] = data
	return table

class dataSet:
	def __init__(self, set1, set2, count):
		self.bset = set1
		self.allItemNames = set2
		self.count = count
	def __hash__(self):
		return hash(frozenset(self.bset))
	def return_vars(self):
		return self.bset, self.allItemNames
	def return_count(self):
		return self.count
	def __eq__(self, other):
		return (self.__class__ == other.__class__ and self.bset == other.bset)

class rule:
	def __init__(self, head, body):
		self.head = head
		self.body = body
	def __hash__(self):
		return hash(frozenset(self.head))+hash(frozenset(self.body))
	def return_vars(self):
		return self.body, self.head
	def return_rule(self):
		return self.head.union(self.body)
	def return_all(self):
		return self.return_rule(), self.body, self.head
	def __eq__(self, other):
		return (self.__class__ == other.__class__ and self.head == other.head and self.body == other.body)
	def __str__(self):
		return "%s -> %s\n" % (str(set(map(lambda x: x.decode("utf-8"), self.return_vars()[0]))), str(set(map(lambda x: x.decode("utf-8"), self.return_vars()[1]))))

def freqSetSetup(saveToFile):
	my_data = genfromtxt(testFile, delimiter='\t', dtype="|S10")

	formatted_data = formatData(my_data)
	allItemNames = set()
	np.vectorize(lambda x: allItemNames.add(x))(formatted_data)
	passin = set()
	passin.add(dataSet(set(), allItemNames, len(formatted_data)))
	set0 = aporoi_freq_set(passin, formatted_data, support)

	if(saveToFile):
		with open(freqSetFileName, 'wb') as picfile:
			pickle.dump(set0, picfile, protocol=pickle.HIGHEST_PROTOCOL)
	return set0
	
def rulesSetup(freqSet,saveToFile):

	freqTable = formatFreqSet(freqSet)
	rulesSet = rule_gen(freqTable, confidence)

	for rule in rulesSet:
		print("%s -> %s" % rule.return_vars())

	if(saveToFile):
		with open(ruleFileName, 'wb') as picfile:
			pickle.dump(rulesSet, picfile, protocol=pickle.HIGHEST_PROTOCOL)
	return rulesSet


class AssoRule:
	def __init__(self, ruleSet):
		self.ruleSet = ruleSet
	def template1(self, ide, quantity, array):
		index = 0
		pred = None
		predCount = 0
		if(ide == "BODY"):
			index = 1
		elif(ide == "HEAD"):
			index = 2
		else:
			assert ide == "RULE"
		if(quantity == "ANY"):
			#print(len(self.ruleSet))
			#list(map(lambda rulez: print(rulez.return_all()[index], array), self.ruleSet))
			pred = lambda rulez: len(set(map(lambda x: x.decode("utf-8"),rulez.return_all()[index])).intersection(set(array))) > 0
		elif(quantity != "NONE"):
			predCount = quantity
		if(pred == None):
			pred = lambda rulez: len(set(map(lambda x: x.decode("utf-8"),rulez.return_all()[index])).intersection(set(array))) == quantity
		resultSet = set(filter( pred, self.ruleSet))
		return resultSet, len(resultSet)

	def template2(self, ide, quantity):
		index = 0
		if(ide == "BODY"):
			index = 1
		elif(ide == "HEAD"):
			index = 2
		else:
			assert ide == "RULE"
		resultSet = set(filter( lambda rulez: len(rulez.return_all()[index]) == quantity, self.ruleSet))
		return resultSet, len(resultSet)
	def template3(self, ide, *arbituraryData):
		t0, t1= None, None
		union = False
		if("or" in ide):
			[t0, t1] = ide.split("or");
			union = True
		else:
			assert "and" in ide
			[t0, t1] = ide.split("and");
		set0, set1, resultSet = None, None, None

		assert (t0== "1" or t0 == "2")
		set0, cnt =  self.template1(arbituraryData[0], arbituraryData[1], arbituraryData[2]) if (t0 == "1") else self.template2(arbituraryData[0], arbituraryData[1])
		arbituraryData = arbituraryData[3:] if (t0 == "1") else arbituraryData[2:]
		
		assert (t1== "1" or t1 == "2")
		set1, cnt = self.template1(arbituraryData[0], arbituraryData[1], arbituraryData[2]) if (t1 == "1") else self.template2(arbituraryData[0], arbituraryData[1])
		
		if(union):
			resultSet = set0.union(set1)
		else:
			resultSet = set0.intersection(set1)
		return resultSet, len(resultSet)
			
# start of main program
freqSet = None
try:
	with open(freqSetFileName, 'rb') as picfile:
		freqSet = pickle.load(picfile)
except FileNotFoundError:
	freqSet = freqSetSetup(saveToFile = True)

with open(freqPrintFileName, 'w') as picfile:
	#for rulez in freqSet:	
	#	picfile.write("%s\n" % str(set(map(lambda x: x.decode("utf-8"), rulez.return_vars()[0]))))
	i = 0
	picfile.write("Support is set to be: %f\n" % (support*100))
	while True:
		i+=1
		list0 = list(filter(lambda x: len(x.return_vars()[0]) == i,freqSet))
		if(len(list0) == 0):
			picfile.write("number of all lengths frequent itemsets: %d\n" % len(freqSet))
			break;
		else:
			picfile.write("number of length-%d lengths frequent itemsets: %d\n" % (i, len(list0)))

freqRule = None
try:
	with open(ruleFileName, 'rb') as picfile:
		freqRule = pickle.load(picfile)
except FileNotFoundError:
	freqRule = rulesSetup(freqSet, saveToFile = True)

with open(rulePrintFileName, 'w') as ruleFile:
	for rule in freqRule:
		ruleFile.write("%s -> %s\n" % (str(set(map(lambda x: x.decode("utf-8"), rule.return_vars()[0]))), str(set(map(lambda x: x.decode("utf-8"), rule.return_vars()[1])))))

asso_rule = AssoRule(freqRule)

result11, cnt11 = asso_rule.template1("RULE", "ANY", ['G59_Up'])
result12, cnt12 = asso_rule.template1("RULE", "NONE", ['G59_Up'])
result13, cnt13 = asso_rule.template1("RULE", 1, ['G59_Up'])
result14, cnt14 = asso_rule.template1("HEAD", "ANY", ['G59_Up'])
(result15, cnt15) = asso_rule.template1("HEAD", "NONE", ['G59_Up'])
(result16, cnt16) = asso_rule.template1("HEAD", 1, ['G59_Up', 'G10_Down'])
(result17, cnt17) = asso_rule.template1("BODY", "ANY", ['G59_Up'])
(result18, cnt18) = asso_rule.template1("BODY", "NONE", ['G59_Up'])
(result19, cnt19) = asso_rule.template1("BODY", 1, ['G59_Up', 'G10_Down'])

(result21, cnt21) = asso_rule.template2("RULE", 3)
(result22, cnt22) = asso_rule.template2("HEAD", 2)
(result23, cnt23) = asso_rule.template2("BODY", 1)

(result31, cnt31) = asso_rule.template3("1or1", "HEAD", "ANY", ['G10_Down'],
"BODY", 1, ['G59_Up'])
(result32, cnt32) = asso_rule.template3("1and1", "HEAD", "ANY",
['G10_Down'], "BODY", 1, ['G59_Up'])
(result33, cnt33) = asso_rule.template3("1or2", "HEAD", "ANY", ['G10_Down'],
"BODY", 2)
(result34, cnt34) = asso_rule.template3("1and2", "HEAD", "ANY",
['G10_Down'], "BODY", 2)
(result35, cnt35) = asso_rule.template3("2or2", "HEAD", 1, "BODY", 2)
(result36, cnt36) = asso_rule.template3("2and2", "HEAD", 1, "BODY", 2)

finalList = ["result11", "result12", "result13", "result14", "result15", "result16", "result17", "result18", "result19", "result21", "result22", "result23", "result31", "result32", "result33", "result34", "result35", "result36"]

finalcnt = [cnt11, cnt12, cnt13, cnt14, cnt15, cnt16, cnt17, cnt18, cnt19, cnt21, cnt22, cnt23, cnt31, cnt32, cnt33, cnt34, cnt35, cnt36]

with open(ruleTemplatePrintFileName, 'w') as tempFile:
	for index, listthing in enumerate(finalList):
		tempFile.write(listthing + ":\n")
		tempFile.write("count: %d\n"%(finalcnt[index]))
		for rulez in eval(listthing):
			tempFile.write(str(rulez))


#print(aporoi_freq_set(formatted_data, [b"1 Up", b"2 Down"], 0.6))
