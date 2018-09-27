import numpy as np
from collections import deque


class Node:
    def __init__(self, cols, sample, rows):
        self.cols = cols
        self.last_idx = cols[-1] + 1
        self.avaiable_rows = rows
        self.sample = sample

    def __repr__(self):
        return "idx: {} \n ele: {}\n".format(str(self.cols), str(self.sample))

    def __str__(self):
        return "idx: {} \n ele: {}\n".format(str(self.cols), str(self.sample))


class apriori:
    def __init__(self, filename):
        self.data = np.genfromtxt(filename, delimiter='\t', dtype=str)
        self.freqset = {}
        self.gen = {}
        self.num_gen = len(self.data[0]) - 1
        self.freqlen = np.zeros(self.num_gen, dtype=int)

    def apriori(self, sup):
        self.reset()
        self.init_first_level(sup)
        i = 0
        while len(self.freqset[i]) > 0:
            temp_node = self.freqset[i].popleft()
            for n in range(temp_node.last_idx, self.num_gen, 1):
                cols = self.obtain_gen_seq([n])
                up = []
                down = []
                for elem in temp_node.avaiable_rows:
                    if cols[elem] == "Up":
                        up.append(elem)
                    else:
                        down.append(elem)

                if len(up) / self.num_gen >= sup:
                    new_cols = temp_node.cols.copy()
                    new_cols.append(n)
                    self.freqset[i + 1].append(Node(new_cols, cols[up[0]], up))
                if len(down) / self.num_gen >= sup:
                    new_cols = temp_node.cols.copy()
                    new_cols.append(n)
                    self.freqset[i + 1].append(Node(new_cols, cols[down[0]], down))

            if len(self.freqset[i]) == 0 and len(self.freqset[i + 1]) != 0:
                i += 1
                self.freqlen[i] = len(self.freqset[i])

        print(self.freqlen)

    def init_first_level(self, sup):
        for i in range(self.num_gen):
            cols = self.obtain_gen_seq([i])
            self.freqset[i] = deque()
            up = []
            down = []
            for idx in range(len(cols)):
                if cols[idx] == "Up":
                    up.append(idx)
                else:
                    down.append(idx)
            if len(up) / self.num_gen >= sup:
                self.freqset[0].append(Node([i], cols[up[0]], up))
            if len(down) / self.num_gen >= sup:
                self.freqset[0].append(Node([i], cols[down[0]], down))

        self.freqlen[0] = len(self.freqset[0])

    def obtain_gen_seq(self, seq_col_idx, seq_row_idx=None):
        if seq_row_idx is None:
            return self.data[:, seq_col_idx]

        return self.data[seq_row_idx, seq_col_idx]

    def reset(self):
        self.freqset = {}
        self.gen = {}
        self.num_gen = len(self.data[0]) - 1
        self.freqlen = np.zeros(self.num_gen, dtype=int)

    def output_ans(self, sup):
        with open("freqset{}.txt".format(sup), "w") as f:
            f.write("Support is set to be {}%\n".format(sup * 100))
            for idx in range(len(self.freqlen)):
                f.write("number of length-{} frequent itemsets: {}\n".format(idx + 1, self.freqlen[idx]))
            f.close()


if __name__ == '__main__':
    ap = apriori("associationruletestdata.txt")
    ap.apriori(0.3)
    ap.output_ans(0.3)
