import numpy as np
from numpy import linalg as la
from sympy import *


############### nonlinear equation ###############


x = symbols('x')

iterative_formula = 1 / (2 ** x)
iterative_formula_diff = diff(iterative_formula)

IterativeFuntion = lambdify(x, iterative_formula)
IterativeFuntionDiff = lambdify(x, iterative_formula_diff)

ksi = 0.5 
Q = IterativeFuntionDiff(ksi) # to evaluate the condition of achieving accuracy


def AccuracyCondition(x_n_plus1, x_n, q = Q):
    return abs(x_n_plus1 - x_n) / (1 - q)


def SimpleIterativeMethod(x_0, epsilon = 1e-16):
    res = x_0
    tmp = res - 1 # -1 because there is no do-while in python((
    iters = 0
    while AccuracyCondition(res, tmp) > epsilon:
        tmp = res
        res = IterativeFuntion(tmp)
        iters += 1
    return res, iters


############### system of nonlinear equations ###############


def F(x):
    F = np.array([[0.0], [0.0]])
    F[0] = np.sin(x[0] + 1.0) - x[1] - 1.2
    F[1] = 2 * x[0] + np.cos(x[1]) - 2.0
    return F

def J(x):
    J = np.identity(2)
    J[0, 0] = np.cos(x[0] + 1.0)
    J[0, 1] = -1.0
    J[1, 0] = 2.0
    J[1, 1] = -np.sin(x[1])
    return J

def Norm3(x):
    return la.norm(x, 2)

def NewtonMethod(startPoint, F, J, epsilon, iterations):
    res = startPoint

    for i in range(0, iterations):
        temp   = res
        res = res - np.matmul(np.linalg.inv(J(res)), F(res))
        if(Norm3(res - temp) < epsilon):
            return res, i

    return res


############### main ###############


def main():
    res, iters = SimpleIterativeMethod(1, 1e-14)
    print("nonlinear equation:", "x =", res, ";", "number of iteraitions:", iters, ";", "epsilon:", epsilon)
    res, iters = NewtonMethod(np.array([[0], [0]]), F, J, 1e-3, 500)
    print("system of nonlinear equations:", "x =", res, ";", "number of iteraitions:", iters, ";", "epsilon:", 1e-3)
    return

if __name__ == '__main__':
    main()