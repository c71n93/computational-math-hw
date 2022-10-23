import numpy as np
from numpy import linalg as la
from itertools import combinations
import scipy


############### task ###############


def FillMatrix():
    n = 10
    A = np.identity(n)
    for i in range(0, n):
        for j in range(0, n):
            if i == j:
                A[i][j] = 1
            else:
                A[i][j] = 1 / (i + j + 2)
    return A

def FillVectorF():
    n = 10
    f = np.array([[1 / i] for i in range(1, n + 1)])
    return f


############### direct-solution-method ###############


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


def GetLU(A):
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i <= j:
                U[i, j] = A[i, j]
                for k in range(i):
                    U[i, j] -= L[i, k] * U[k, j]
            else:
                L[i, j] = 1 / U[j, j] * A[i, j]
                for k in range(j):
                    L[i, j] -= 1 / U[j, j] * L[i, k] * U[k, j]
        
    return L, U


def SolveUpperTriangle(UT, b):
    n = UT.shape[0]
    x = np.zeros(n)

    x[n - 1] = b[n - 1] / UT[n - 1][n - 1]

    for k in range(n - 2, -1, -1):
        sum = 0
        for i in range(k + 1, n):
            sum += UT[k][i] * x[i]
        x[k] = (b[k] - sum) / UT[k][k]

    return x.reshape([n, 1])


def SolveLowerTriangle(LT, b):
    n = LT.shape[0]
    x = np.zeros(n)

    x[0] = b[0] / LT[0][0]

    for k in range(1, n):
        sum = 0
        for i in range(0, k):
            sum += LT[k][i] * x[i]
        x[k] = (b[k] - sum) / LT[k][k]

    return x.reshape([n, 1])


def SolveLUDecomposition(A, b):
    LU = GetLU(A)
    y = SolveLowerTriangle(LU[0], b)
    x = SolveUpperTriangle(LU[1], y)
    residual = Norm3(np.matmul(A, x) - b)
    print("LU decomposition residual:", residual, end="\n\n")
    return x


############### iteration-solution-method ###############


def SolveUpperRelaxation(A, f, omega = 1.5, convergence_criteria = 1e-10):
    init_vector=np.array([[0] for i in range(0, A.shape[0])])
    n = A.shape[0]
    
    D = np.zeros([n, n])
    L = np.zeros([n, n])
    U = np.zeros([n, n])

    x = init_vector

    for i in range(0, n):
        for j in range(0, n):
            if i > j:
                L[i, j] = A[i, j]
            elif i < j:
                U[i, j] = A[i, j]
            else:
                D[i, j] = A[i, j]

    B = np.linalg.inv(D + omega * L)
    C = np.matmul(B, (omega - 1) * D + omega * U)

    while Norm3(f - np.matmul(A, x)) >= convergence_criteria:
        x = -np.matmul(C, x) + np.matmul(omega * B, f)
    
    residual = Norm3(np.matmul(A, x) - f)
    print("Upper Relaxation convergence criterion:", convergence_criteria)
    print("Upper Relaxation residual:", residual, end="\n\n")

    return (x)


############### norms-and-other-functions ###############


def Normalize(x):
    norm = Norm3(x)
    if norm == 0:
        return np.array([[0] for i in range(0, x.shape[0])])
    else:
        return x / norm


def MaxEigenvalue(A):
    iterations = 500
    maxEigenValue = 0.0
    prevVector = np.array([[i + 1] for i in range(0, A.shape[0])])
    currentVector = prevVector

    for i in range(0, iterations):
        prevVector = currentVector
        currentVector = np.matmul(A, currentVector)

    maxEigenValue = np.max(np.abs(currentVector[5] / prevVector[5]))

    return maxEigenValue


def MinEigenvalue(A):
    return 1 / MaxEigenvalue(la.inv(A))


def Norm1(x):
    return la.norm(x, np.inf)


def Norm2(x):
    return la.norm(x, 1)


def Norm3(x):
    return la.norm(x, 2)


def PrintConditionNumbers(A):
    A_inv = la.inv(A)
    m1 = Norm1(A) * Norm1(A_inv)
    m2 = Norm2(A) * Norm2(A_inv)
    m3 = Norm3(A) * Norm3(A_inv)
    print("Condition number in norm 1:", m1)
    print("Condition number in norm 2:", m2)
    print("Condition number in norm 3:", m3, end="\n\n")


def PrintMinMaxEigenvalues(A):
    print("Max eigenvalue:", MaxEigenvalue(A))
    print("Min eigenvalue:", MinEigenvalue(A), end="\n\n")


############### tests ###############


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


############### main ###############


def main():

    A = FillMatrix()
    f = FillVectorF()
    PrintConditionNumbers(A)
    PrintMinMaxEigenvalues(A)
    SolveLUDecomposition(A, f)
    SolveUpperRelaxation(A, f)
    

if __name__ == '__main__':
    main()