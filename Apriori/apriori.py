import numpy as np
from collections import deque
from timeit import default_timer as timer


class Head:
    def __init__(self, elem=[]):
        self.elem = elem

    def __len__(self):
        return len(self.elem)


class Body:
    def __init__(self, elem=[]):
        self.elem = elem

    def __len__(self):
        return len(self.elem)


class Rule:
    def __init__(self, head=None, body=None):
        self.head = head
        self.body = body
        self.elem = head.elem.copy()
        self.elem.extends(body.elem)

    def __len__(self):
        return len(self.head) + len(self.head)


class Node:
    def __init__(self, parent, cols, sample, rows):
        self.cols = cols
        self.last_idx = cols[-1] + 1
        self.avaiable_rows = rows
        if parent is not None:
            self.sample = parent.sample.copy()
            self.sample.extend(sample)
        else:
            self.sample = sample
            # print(self)

    def __repr__(self):
        return "idx: {} \n ele: {}\n".format(str(self.cols), str(self.sample))

    def __str__(self):
        return "idx: {} \n ele: {}\n".format(str(self.cols), str(self.sample))


class apriori:
    def __init__(self, filename):
        self.data = np.genfromtxt(filename, delimiter='\t', dtype=str)
        self.freqset = {}
        self.num_patient = len(self.data)
        self.num_gen = len(self.data[0])
        print("number of patient: {}\n number of gen: {}".format(self.num_patient, self.num_gen))
        self.freqlen = np.zeros(self.num_gen, dtype=int)

    def apriori(self, sup, conf):
        self.reset()
        self.gen_freq(sup)
        self.gen_rule(conf)

    def gen_freq(self, sup):
        self.reset()
        self.init_first_level(sup)
        i = 0
        j = 0  # index of i's length freq set
        while len(self.freqset[i]) > 0:

            temp_node = self.freqset[i][j]
            for n in range(temp_node.last_idx, self.num_gen, 1):
                cols = self.obtain_gen_seq([n])  # time: 6 mins

                args = {}
                for elem in temp_node.avaiable_rows:
                    col = str(cols[elem][0])

                    if col not in args:

                        args[col] = [elem]

                    else:
                        args[col].append(elem)

                for arg in args:

                    if len(args[arg]) / self.num_patient >= sup:
                        new_cols = temp_node.cols.copy()
                        new_cols.append(n)
                        self.freqset[i + 1].append(Node(temp_node, new_cols, [arg], args[arg]))
            j += 1
            if j == len(self.freqset[i]):
                i += 1
                j = 0
                self.freqlen[i] = len(self.freqset[i])
                if len(self.freqset[i]) == 0 or i == len(self.freqset):
                    break

        self.output_ans(sup)
        print(self.freqlen)

    def gen_rule(self, conf):
        pass

    def init_first_level(self, sup):
        for i in range(self.num_gen):
            cols = self.obtain_gen_seq([i])
            self.freqset[i] = deque()
            args = {}
            for elem in range(len(cols)):
                col = str(cols[elem][0])
                if col not in args:
                    args[col] = [elem]
                else:
                    args[col].append(elem)
            for arg in args:
                # print("arg: {} {}".format(arg,len(args[arg])))
                if len(args[arg]) / self.num_patient >= sup:
                    self.freqset[0].append(Node(None, [i], [arg], args[arg]))
        self.freqlen[0] = len(self.freqset[0])

    def obtain_gen_seq(self, seq_col_idx, seq_row_idx=None):
        if seq_row_idx is None:
            return self.data[:, seq_col_idx]

        return self.data[seq_row_idx, seq_col_idx]

    def reset(self):
        self.freqset = {}
        self.gen = {}
        self.num_gen = len(self.data[0])
        self.freqlen = np.zeros(self.num_gen, dtype=int)

    def output_ans(self, sup):
        with open("freqset{}.txt".format(sup), "w") as f:
            f.write("Support is set to be {}%\n".format(sup * 100))
            for idx in range(len(self.freqlen)):
                f.write("number of length-{} frequent itemsets: {}\n".format(idx + 1, self.freqlen[idx]))
            f.close()


if __name__ == '__main__':
    ap = apriori("associationruletestdata.txt")
    start = timer()
    ap.apriori(0.3, 0.7)
    end = timer()
    print("time elapse: {}".format(end - start))
