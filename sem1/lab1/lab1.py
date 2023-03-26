import numpy as np
import matplotlib.pyplot as plt
from sympy import *


def diff_approx_1 (func, x, h):
    f_derrivate = (func(x + h) - func(x)) / h
    return f_derrivate


def diff_approx_1_m (func, x, h):
    f_derivate = (func(x) - func(x - h)) / h
    return f_derivate


def diff_approx_2 (func, x, h):
    f_derivate = (func(x + h) - func(x - h)) / (2 * h)
    return f_derivate


def diff_approx_3 (func, x, h):
    frac1 = (func(x + h) - func(x - h)) / (2 * h)
    frac2 = (func(x + 2 *  h) - func(x - 2 * h)) / (4 * h)
    f_derivate = (4 / 3) * frac1 - (1 / 3) * frac2
    return f_derivate


def diff_approx_4 (func, x, h):
    frac1 = (func(x + h) - func(x - h)) / (2 * h)
    frac2 = (func(x + 2 *  h) - func(x - 2 * h)) / (4 * h)
    frac3 = (func(x + 3 *  h) - func(x - 3 * h)) / (6 * h)
    f_derivate = (3 / 2) * frac1 - (3 / 5) * frac2 + (1 / 10) * frac3
    return f_derivate


diff_approx = [diff_approx_1, diff_approx_1_m, diff_approx_2, diff_approx_3, diff_approx_4]


def get_delta_array(sym_func, x_0, h_arr, diff_approx_func):
    n_points = len(h_arr)
    delta_arr = [0.0] * n_points
    for i in range(0, n_points):
        x = symbols('x')
        func_lambda = lambdify(x, sym_func)
        func_derivate_lambda = lambdify(x, diff(sym_func))
        aprox_derivate = diff_approx_func(func_lambda, x_0, h_arr[i])
        delta_arr[i] = abs(func_derivate_lambda(x_0) - aprox_derivate)
    return delta_arr


def generate_plot(h_arr, delta_arr, func, x_0, approx_num, plot):
    plot.loglog(h_arr, delta_arr, label="№%d" % approx_num, base=2)
    plot.set_xlabel("$h$")
    plot.set_ylabel("$Delta$")
    plot.set_title("Погрешность для производной " + str(func) + " в точке " + "%d" % x_0)
    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plot.grid()


x = symbols('x')

functions = [sin(x * x), cos(sin(x)), exp(cos(sin(x))), ln(x + 3), (x + 3) ** (1 / 2)]


def main():
    n_points = 21
    h_arr = 2 / 2 ** np.arange(1, n_points + 1)

    x_0 = 100

    for f in functions:
        plt.figure(figsize=(9, 5))
        plot = plt.subplot()
        for i in range(len(diff_approx)):
            delta_arr = get_delta_array(f, x_0, h_arr, diff_approx[i])
            generate_plot(h_arr, delta_arr, f, x_0, i + 1, plot)
        plt.show()


if __name__ == "__main__":
    main()