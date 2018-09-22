import numpy as np


class PCA:
    def __init__(self,*filename):
        data_list = []
        for file in filename:
            data_list.append(self.read_file(file))


    def pca(self):
        pass

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

    # compute mean value of a list of number

    def mean(self, vector: list):
        return sum(vector) / len(vector)

        # calculate variant list

    def compute_variant(self, vector: list, mean: float):
        return [i - mean for i in vector]

    def compute_covariant(self, matrix_x, matrix_covx):
        COV = np.ndarray(shape=(len(matrix_x), len(matrix_covx)), dtype=np.float32)
        for idx in range(len(matrix_x)):
            for co_idx in range(len(matrix_covx)):
                COV[idx][co_idx] = 0
                # print("*****************************************")
                for x, x1 in zip(matrix_x[idx], matrix_covx[co_idx]):
                    # print("x1:{} x2:{}".format(x, x1))
                    # print(x * x1)
                    COV[idx][co_idx] += (x * x1)

                COV[idx][co_idx] = (COV[idx][co_idx] / (len(matrix_x[0])))
        return COV


if __name__ == '__main__':
    pca = PCA()
    matrix_x = [[19, 39, 30, 30, 15, 15, 15, 30], [63, 74, 87, 23, 35, 43, 32, 73]]
    matrix_covx = [[-5.1, 14.9, 5.9, 5.9, -9.1, - 9.1, - 9.1, 5.9],
                   [9.25, 20.25, 33.25, -30.75, -18.75, -10.75, -21.75, 19.25]]
    print(pca.compute_covariant(matrix_x, matrix_covx))
