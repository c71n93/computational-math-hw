import numpy as np
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


############### main ###############


def main():
    x_0 = 1
    epsilon = 1e-14
    res, iters = SimpleIterativeMethod(x_0, epsilon)
    print("x =", res, ";", "number of iteraitions:", iters, ";", "epsilon:", epsilon)

if __name__ == '__main__':
    main()