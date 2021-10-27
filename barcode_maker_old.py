#!/usr/bin python
_author_ = "Mitch Ledwith"
"""
adding this because it's a fun little problem
lets make 24 unique 5-nt barcodes, out of 1024 options
all requiring at least 2 mutations to get to another barcode (i.e. one mutation won't do it)
all between 25-75% GC ???? (no homopolymers) (avoid 3-identical bases in a row)
"""


from itertools import product
from itertools import combinations
from collections import Counter
import random

def onebpdifference(seq): #returns list of barcodes that have one bp difference relative to starting barcode
	base_list = ["A","T","C","G"]
	temp_list = []
	for x in range(0,len(seq)): #swaps each position out for a different nucleotide
		for item in base_list:
			my_new_creation = seq[:x] + item + seq[x+1:]
			if my_new_creation != seq:
				temp_list.append(my_new_creation)
	return temp_list 
			
bases = "ATCG"
permutations = ["".join(i) for i in product(bases,repeat=5)] #generates all possible 5-nt barcodes out of 4 bases

new_perm_list = []

for item in permutations:
	if "AAA" not in item and "GGG" not in item and "CCC" not in item and "TTT" not in item: #gets rid of homopolymers
		my_count = Counter(item)
		if 1 < (float(my_count["A"])+float(my_count["T"])) < 4: #gets rid of anything that is < 40% GC or >60% GC
			new_perm_list.append(item)


random.Random(6).shuffle(new_perm_list) #needs to be shuffled, otherwise favors the first members of the list
one_base_library = []
random.seed() #reseed random

for item in new_perm_list: #filters out barcodes that don't pass the one-bp difference rule after starting with random barcode
	new_list = onebpdifference(item)
	if set(new_list).isdisjoint(set(one_base_library)) == True:
		one_base_library.append(item)



A_list = [i for i in one_base_library if i[0] == "A"] #barcodes with first base == A
T_list = [i for i in one_base_library if i[0] == "T"]
C_list = [i for i in one_base_library if i[0] == "C"]
G_list = [i for i in one_base_library if i[0] == "G"]




def balanced_question_mark(list_codes,x): #check to see if the barcodes at position x have equal base frequency
	new_list = [i[x] for i in list_codes]
	balance = Counter(new_list)
	if balance["A"] == 6 and balance["T"] == 6 and balance["G"] == 6 and balance["C"] == 6:
		return True
	else:
		return False
	

z = 0
Acomb = combinations(A_list,6) #generator function of all possible combinations 6 barcodes that start with A
Tcomb = combinations(T_list,6)
Ccomb = combinations(C_list,6)
Gcomb = combinations(G_list,6)
all_comb = product(Acomb,Tcomb,Ccomb,Gcomb) #generator function of generator function, to avoid storing all of this in memory
for i in all_comb:
	i_not = [j for k in i for j in k]
	balance = balanced_question_mark(i_not,1)
	if balance == True: #if base balance at position 1 is good, check position 2
		balance = balanced_question_mark(i_not,2)
		if balance == True:
			balance = balanced_question_mark(i_not,3)
			if balance == True:
				balance = balanced_question_mark(i_not,4)
				if balance == True:
					print i_not
					z += 1
					print z

		



	
