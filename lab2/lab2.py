import numpy as np
from numpy import linalg as la
from itertools import combinations


def GetLU(A):
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros((n, n))

    for i in range(0, n):
        for j in range(0, n):
            if i <= j:
                sum = 0
                for k in range(0, i):
                    sum += L[i][k] * U[k][j]
                U[i][j] = A[i][j] - sum

            elif i > j:
                sum = 0
                for k in range(0, j):
                    L[i][k] * U[k][j]
                L[i][j] = (A[i][j] - sum) / U[j][j]

    return (L, U)


def SolveUpperTriangle(UT, b):
    n = UT.shape[0]
    x = np.zeros(n)

    x[n - 1] = b[n - 1] / UT[n - 1][n - 1]

    for k in range(n - 2, -1, -1):
        sum = 0
        for i in range(k + 1, n):
            sum += UT[k][i] * x[i]
        x[k] = (b[k] - sum) / UT[k][k]

    return x


def SolveLowerTriangle(LT, b):
    n = LT.shape[0]
    x = np.zeros(n)

    x[0] = b[0] / LT[0][0]

    for k in range(1, n):
        sum = 0
        for i in range(0, k):
            sum += LT[k][i] * x[i]
        x[k] = (b[k] - sum) / LT[k][k]

    return x


def SolveLUDecomposition(A, b):
    LU = GetLU(A)
    y = SolveLowerTriangle(LU[0], b)
    x = SolveUpperTriangle(LU[1], y)
    return x


def IsLUDecomposable(A):
    n = A.shape[0]
    index_arr = np.arange(0, n)

    for order in range(1, n + 1):
        perms = list(combinations(index_arr, order))
        for perm in perms:
            M = A[[[j] for j in perm ], [j for j in perm]]
            if la.det(M) == 0:
                return False

    return True


###############################################################################


def normalize(x): #refactor
    norm = Norm3(x)
    if norm == 0:
        return np.array([[0] for i in range(0, x.shape[0])])
    else:
        return x / norm


def MaxEigenvalue(A): #refactor
    iterations  = 500
    prev_vector = prev_vector = np.array([[i + 1] for i in range(0, A.shape[0])])
    curr_vector = prev_vector

    for i in range(iterations):
        curr_vector = normalize(curr_vector)
        prev_vector = normalize(prev_vector)

        prev_vector = curr_vector
        curr_vector = np.matmul(A, curr_vector)

    max_eigval = np.max(curr_vector / prev_vector)
    return max_eigval


def MinEigenvalue(A): #refactor
    return 1 / MaxEigenvalue(la.inv(A))


def Norm1(x):
    return la.norm(x, np.inf)


def Norm2(x):
    return la.norm(x, 1)


def Norm3(x):
    return la.norm(x, 2)


###############################################################################


def SimpleLUTest():
    A = np.array([[4, 3], [6, 3]])
    print(IsLUDecomposable(A))
    LU = GetLU(A)
    print(SolveLowerTriangle(LU[0], np.array([1, 1])) == np.array([1, -0.5]))
    print(SolveUpperTriangle(LU[1], np.array([1, 1.5]))== np.array([1, -1]))


def SolveLUTest():
    A = np.array([[1, 2, 3], [4, 5, 6], [7, 5, 9]])
    b = np.array([10, 11, 12])
    print(IsLUDecomposable(A))
    print(abs(SolveLUDecomposition(A, b) - np.array([-4.5, 0, 29 / 6])) < 10e-15)


def EigenvaluesTest():
    A = np.array([[1, 2, 3], [4, 5, 6], [7, 5, 9]])
    print(MaxEigenvalue(A))
    print(MinEigenvalue(A))
    print(la.eig(A)[0])


###############################################################################


def main():
    SimpleLUTest()
    SolveLUTest()
    EigenvaluesTest()
    

if __name__ == '__main__':
    main()