import numpy as np
import matplotlib.pyplot as plt

B_GLOBAL = 1.5

# (y_1', y_2') = f(y_1, y_2)
def Brusselator(t, y, A = 1.0, B = B_GLOBAL):
  # print("Brusselator:", y)
  f = np.array([1.0, 0.0])
  f[0] = A + y[0] * y[0] * y[1] - (B + 1) * y[0]
  f[1] = B * y[0] - y[0] * y[0] * y[1]
  # print(f[0], f[1])
  return f

# y' = f(t, y)
def RungeKuttaIteration_4Order(t_n, y_n, h, f):
  # print("RungeKuttaIteration_4Order:", y_n)
  f_1 = f(t_n, y_n)
  f_2 = f(t_n + h / 2, y_n + f_1 * h / 2)
  f_3 = f(t_n + h / 2, y_n + f_2 * h / 2)
  f_4 = f(t_n + h, y_n + f_3 * h)
  y_n_1 = y_n + h / 6 * (f_1 + 2 * f_2 + 2 * f_3 + f_4)
  # print("RungeKuttaIteration_4Order (y_n+1):", y_n_1)
  return y_n_1

def BrusselatorSolutionsRK(N = int(1e4), h = 0.01):
  start = np.array([1.0, 1.0])
  t_array = [0.0]
  y_array = [start]

  for i in range(0, N):
    # print("\n\nBrusselatorSolutionsRK:", y_array[i])
    y_new = RungeKuttaIteration_4Order(t_array[i], y_array[i], h, Brusselator)
    # print("BrusselatorSolutionsRK (y_new):", y_new)
    y_array.append(y_new)
    # print("BrusselatorSolutionsRK (y_array):", y_array)
    t_array.append(h * (i + 1))

  return y_array, t_array

def ShowPlotsRK(u, v, t):
  plt.figure(figsize=[12, 5], dpi=100)
  plt.plot(u, v)
  plt.title("Рунге-Кутта. Фазовый портрет v(u) при B = %.2f" %B_GLOBAL)
  plt.grid()
  plt.show()

  plt.figure(figsize=[12, 5], dpi=100)
  plt.plot(t, u, label = "u")
  plt.plot(t, v, label = "v")
  plt.xlabel("t")
  plt.ylabel("u, v")
  plt.legend()
  plt.title("Рунге-Кутта. u(t) и v(t) при B = %.2f" %B_GLOBAL)
  plt.grid()
  plt.show()

def main():
  # print(Brusselator(t = 1, y = [0.9975, 1.0025], B = 1.5))
  N = 10000
  RK_solution_arr, t_arr = BrusselatorSolutionsRK(N = N)
  RK_solution_arr_transposed = np.array(RK_solution_arr).transpose()
  ShowPlotsRK(RK_solution_arr_transposed[0], RK_solution_arr_transposed[1], t_arr)

if __name__ == "__main__":
    main()