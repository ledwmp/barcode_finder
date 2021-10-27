from itertools import product
from itertools import combinations

class barcode:
    """Class to store barcode and generate barcode stats

    """
    def __init__(self,seq):
        self.seq = seq
        self.gc = (seq.count("A")+seq.count("T"))/float(len(seq))
    def one_bp_diff(self):
        """Method to generate list of barcodes with one mutation away from input barcode

        """
        base_list = ["A","T","C","G"]
        diff_list = []
        for i in range(0,len(self.seq)):
            for base in base_list:
                altbc = self.seq[:i] + base + self.seq[i+1:]
                if altbc != self.seq:
                    diff_list.append(altbc)
        return diff_list
    def pass_filters(self):
        """Method to check for homopolymers and high or low gc-content

        """
        if 0.25 < self.gc < 0.75 and \
        "A"*3 not in self.seq and \
        "G"*3 not in self.seq and \
        "T"*3 not in self.seq and \
        "C"*3 not in self.seq:
            return True
        else:
            return False

class barcode_population:
    """Class to generate population of barcodes, filter members by barcode stats, subset a group of barcodes that have have hamming distance > 1, and then iterate through combinations of these to find a balanced set.

    """
    def __init__(self,membern,bclen,nucs):
        self._membern = membern
        self._bclen = bclen
        self._solutions = []
        self._nucs = nucs
        self._members = []
    def generate_members(self):
        self._members = [barcode("".join(i)) for i in product(self._nucs,repeat=self._bclen)]
        return self._members
    def filter_members(self):
        self._filtermembers = [i for i in self.generate_members() if i.pass_filters() == True]
        return self._filtermembers
    def hammed_members(self):
        tmp_list = []
        for i in self.filter_members():
            if set(i.one_bp_diff()).isdisjoint(set([j.seq for j in tmp_list])) == True:
                tmp_list.append(i)
        self._hammedmembers = tmp_list
        return self._hammedmembers
    def separate_firstnuc(self):
        keys = [i for i in self._nucs]
        in_list = self.hammed_members()
        out_dict = {i:[j for j in in_list if j.seq[0] == i] for i in keys}
        self._firstbasedict = out_dict
        return self._firstbasedict
    def nuc_combos(self):
        keys = [i for i in self._nucs]
        in_dicts = self.separate_firstnuc()
        combo_dict = {i:combinations(in_dicts[i],self._membern/len(self._nucs)) for i in keys}
        all_combinations = product(*combo_dict.values())
        self._allcombos = all_combinations
        return self._allcombos
    def balance_check(self,tmp_list):
        bc_len = self._bclen
        nuc_list = [i for i in self._nucs]
        tmp_list = [i.seq for j in tmp_list for i in j]
        for i in range(0,bc_len):
            pos_list = "".join([j[i] for j in tmp_list])
            count_list = [pos_list.count(j) for j in nuc_list]
            if count_list.count(count_list[0]) != len(count_list):
                return False
        return True

def main():
    bcpop = barcode_population(24,4,"ATCG")
    z = 0
    for i in bcpop.nuc_combos():
        z+=1
        print z
        if bcpop.balance_check(i) == True:
            print "Found a solution"
            print [j.seq for k in i for j in k]
            break
if __name__ == "__main__":
    main()
