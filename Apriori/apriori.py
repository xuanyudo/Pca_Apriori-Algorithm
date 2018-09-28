import numpy as np
from collections import deque
from timeit import default_timer as timer


class Head:
    def __init__(self, elem=[]):
        self.elem = elem

    def __len__(self):
        return len(self.elem)

    def __repr__(self):
        return "head: {} \t length: {}\n".format(list(map(lambda x: x + 1, self.elem)), len(self))

    def __str__(self):
        return "head: {} \t length: {}\n".format(list(map(lambda x: x + 1, self.elem)), len(self))


class Body:
    def __init__(self, elem=[]):
        self.elem = elem

    def __len__(self):
        return len(self.elem)

    def __repr__(self):
        return "body: {} \t length: {}\n".format(list(map(lambda x: x + 1, self.elem)), len(self))

    def __str__(self):
        return "body: {} \t length: {}\n".format(list(map(lambda x: x + 1, self.elem)), len(self))


class Rule:
    def __init__(self, head=None, body=None, node=None):
        self.head = head
        self.body = body
        self.node = node
        self.elem = head.elem.copy()
        self.elem.extend(body.elem)
        self.rule_map = {}
        self.init_rule_map()

    def __len__(self):
        return len(self.head) + len(self.head)

    def __eq__(self, other):
        return (
            self.__class__ == other.__class__ and self.head.elem == other.head.elem and self.body.elem == other.body.elem and self.node == other.node)

    def __hash__(self):
        return hash(frozenset(self.head.elem)) ^ hash(frozenset(self.head.elem)) ^ hash(self.node)

    def __repr__(self):
        return "rule:{} \t {}\n".format(self.head, self.body)

    def __str__(self):
        return "rule:{} \t {}\n".format(self.head, self.body)

    def init_rule_map(self):
        self.rule_map["HEAD"] = []
        self.rule_map["BODY"] = []
        self.rule_map["RULE"] = []
        for key in self.head.elem:
            self.rule_map["HEAD"].append("G{}_{}".format(key + 1, self.node.sample[self.node.cols.index(key)].upper()))
            self.rule_map["RULE"].append("G{}_{}".format(key + 1, self.node.sample[self.node.cols.index(key)].upper()))
        for key in self.body.elem:
            self.rule_map["BODY"].append("G{}_{}".format(key + 1, self.node.sample[self.node.cols.index(key)].upper()))
            self.rule_map["RULE"].append("G{}_{}".format(key + 1, self.node.sample[self.node.cols.index(key)].upper()))

    def repl_elem(self):
        self.head.elem = self.rule_map["HEAD"]
        self.body.elem = self.rule_map["BODY"]


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

    def __hash__(self):
        return hash(frozenset(self.cols)) ^ hash(self.sample[-1])

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and self.cols == other.cols)


class apriori:
    def __init__(self, filename):
        self.data = np.genfromtxt(filename, delimiter='\t', dtype=str)
        self.freqset = {}
        self.num_patient = len(self.data)
        self.num_gen = len(self.data[0])
        self.node_dict = {}
        self.conf_rules = {}
        print("number of patient: {}\n number of gen: {}".format(self.num_patient, self.num_gen))
        self.freqlen = np.zeros(self.num_gen, dtype=int)

    def apriori(self, sup, conf):
        self.reset()
        self.gen_freq(sup)
        self.gen_rule(conf)
        # for rule in self.conf_rules:
        #     rule.repl_elem()
        #     print(rule)
        print(self.template1("RULE", "ANY", ['G1_DOWN']))
        # print(self.template3('1or1', "RULE", "ANY", ['G59_UP'], "RULE", "NONE", ['G1_DOWN']))

    def template1(self, field, number, gen):
        count = 0
        result = []
        if number == "ANY":
            match_need = 1
        elif number == "NONE":
            match_need = 0
        else:
            match_need = number

        for rule in self.conf_rules:
            match = 0
            for elem in gen:
                if elem in rule.rule_map[field]:
                    if number == "ANY":
                        match += 1
                        break
                    else:
                        match += 1
            if match_need == match:
                result.append(rule)
                count += 1
        return result, count

    def template2(self, field, size):
        count = 0
        result = []
        for rule in self.conf_rules:
            if len(rule.rule_map[field]) == size:
                result.append(rule)
                count += 1
        return result, count

    def template3(self, templates, *args):
        t1 = templates[0]
        t2 = templates[-1]
        c1 = 0
        r1 = []
        c2 = 0
        r2 = []
        op = templates[1:-1]
        arg_idx = 0
        if t1 == '1':
            r1, c1 = self.template1(args[arg_idx], args[arg_idx + 1], args[arg_idx + 2])
            arg_idx += 3
        else:
            r1, c1 = self.template2(args[arg_idx], args[arg_idx + 1])
            arg_idx += 2
        if t2 == '1':
            r2, c2 = self.template1(args[arg_idx], args[arg_idx + 1], args[arg_idx + 2])
            arg_idx += 3
        else:
            r2, c2 = self.template2(args[arg_idx], args[arg_idx + 1])
            arg_idx += 2
        if op == 'or':
            res = set(r1).union(set(r2))
            count = len(res)
        else:
            res = set(r1).intersection(set(r2))
            count = len(res)

        return res, count

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
                        new_node = Node(temp_node, new_cols, [arg], args[arg])
                        self.node_dict[new_node] = new_node
                        self.freqset[i + 1].append(new_node)
            j += 1
            if j == len(self.freqset[i]):
                i += 1
                j = 0
                self.freqlen[i] = len(self.freqset[i])
                if len(self.freqset[i]) == 0 or i == len(self.freqset):
                    break

        self.output_ans(sup)
        print(self.freqlen)

        # self.node_dict[Node(None,[9],['Up'],None)]

    def gen_rule(self, conf):
        queue = deque()
        for key in self.node_dict:
            rule = Rule(Head(key.cols), Body([]), key)
            queue.append(rule)
            print(key)
        while len(queue) != 0:
            temp = queue.popleft()
            if len(temp.head.elem) > 1:
                for i in range(len(temp.head.elem)):
                    new_head = temp.head.elem.copy()
                    body = temp.body.elem.copy()
                    body.append(new_head.pop(i))
                    new_head.sort()
                    body.sort()
                    temp_node = self.node_dict[
                        Node(None, new_head, [temp.node.sample[temp.node.cols.index(new_head[-1])]], None)]
                    if len(temp.node.avaiable_rows) / len(temp_node.avaiable_rows) >= conf:
                        rule_conf = Rule(Head(new_head), Body(body), temp.node)
                        if rule_conf not in self.conf_rules:
                            self.conf_rules[rule_conf] = rule_conf
                            queue.append(rule_conf)
                            # else:
                            #     continue

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
                    new_node = Node(None, [i], [arg], args[arg])
                    self.node_dict[new_node] = new_node
                    self.freqset[0].append(new_node)
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
    ap.apriori(0.5, 0.7)
    end = timer()
    print("time elapse: {}".format(end - start))
