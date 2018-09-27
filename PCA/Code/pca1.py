import numpy as np
import matplotlib.pyplot as plt


class PCA:
    def __init__(self):
        self.data_list = []
        self.colorMap = {}
        self.color = ['bo', 'ro', 'go', 'co', 'mo', 'yo', 'ko']
        self.feature_list = []
        self.disease = []

    def pca(self, filename):
        self.data_list = self.read_file(filename)
        self.data_list = np.array(self.data_list)
        self.parse_feature_disease()
        variant_vec = [self.compute_variant(self.feature_list[i], self.mean(self.feature_list[i])) for i in
                       range(len(self.feature_list))]
        matrix_cov = self.compute_covariant(self.feature_list, variant_vec)
        value, vector = self.eigen(matrix_cov)
        y1 = self.fn(vector, value, variant_vec, 0)
        y2 = self.fn(vector, value, variant_vec, 1)
        self.draw_graph(y1, y2)
        # print("value: {} vector: {}".format(value, vector))

    def read_file(self, filename):
        ret = []
        with open(filename, "r") as f:
            line = f.readline()
            while line != "":
                data = line.rstrip('\n').split('\t')
                ret.append(data)
                line = f.readline()
            f.close()
        return ret

    def parse_feature_disease(self):
        self.disease = []
        self.feature_list = []
        for i in range(len(self.data_list[0])):
            if i == len(self.data_list[0]) - 1:
                self.disease.append(self.data_list[:, i:i + 1].flatten())
            else:
                self.feature_list.append([float(x) for x in self.data_list[:, i:i + 1].flatten()])

    # compute mean value of a list of number

    def mean(self, vector: list):
        return sum(vector) / len(vector)

        # calculate variant list

    def compute_variant(self, vector: list, mean: float):
        return [i - mean for i in vector]

    def compute_covariant(self, matrix_x, matrix_varx):
        COV = np.ndarray(shape=(len(matrix_x), len(matrix_varx)), dtype=np.float32)
        for idx in range(len(matrix_x)):
            for co_idx in range(len(matrix_varx)):
                COV[idx][co_idx] = 0
                # print("*****************************************")
                for x, x1 in zip(matrix_x[idx], matrix_varx[co_idx]):
                    # print("x1:{} x2:{}".format(x, x1))
                    # print(x * x1)
                    COV[idx][co_idx] += (x * x1)

                COV[idx][co_idx] = (COV[idx][co_idx] / (len(matrix_x[0])))
        self.eigen(COV)
        return COV

    def eigen(self, cov):
        return np.linalg.eig(cov)

    def fn(self, eigen_vector, eigen_value: np.ndarray, x, order):
        y = []
        idx = eigen_value.argsort()[::-1][order]
        for i in range(len(x[0])):
            value = 0
            for j in range(0, len(x)):
                value += x[j][i] * eigen_vector[idx][j]
            y.append(value)
        return y

    def draw_graph(self, y1, y2):
        self.colorMap = {}
        for i in range(len(self.disease[0])):
            if self.disease[0][i] in self.colorMap.keys():
                plt.plot(y1[i], y2[i], self.colorMap[self.disease[0][i]])
            else:
                self.colorMap[self.disease[0][i]] = self.color[len(self.colorMap)]
                plt.plot(y1[i], y2[i], self.colorMap[self.disease[0][i]])
        plt.show()


if __name__ == '__main__':
    pca = PCA()
    pca.pca("pca_a.txt")
    pca.pca("pca_b.txt")
    pca.pca("pca_c.txt")
