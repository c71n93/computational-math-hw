import numpy as np
import matplotlib.pyplot as plt

B_GLOBAL = 1

# ---------- Condition ----------

# (y_1', y_2') = f(y_1, y_2)
def Brusselator(t, y, A = 1.0, B = B_GLOBAL):
  f = np.array([1.0, 0.0])
  f[0] = A + y[0] * y[0] * y[1] - (B + 1) * y[0]
  f[1] = B * y[0] - y[0] * y[0] * y[1]
  return f

# ---------- Runge-Kutta method ----------

def RungeKuttaIteration_4Order(t_n, y_n, h, f):
  f_1 = f(t_n, y_n)
  f_2 = f(t_n + h / 2, y_n + f_1 * h / 2)
  f_3 = f(t_n + h / 2, y_n + f_2 * h / 2)
  f_4 = f(t_n + h, y_n + f_3 * h)

  y_n_1 = y_n + h / 6 * (f_1 + 2 * f_2 + 2 * f_3 + f_4)
  return y_n_1

def BrusselatorSolutionsRK(N = int(1e4), h = 0.01):
  start = np.array([1.0, 1.0])
  t_array = [0.0]
  y_array = [start]

  for i in range(0, N):
    y_new = RungeKuttaIteration_4Order(t_array[i], y_array[i], h, Brusselator)
    y_array.append(y_new)
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

# ---------- Adams method ----------

def AdamsIteration_4Order(t_prev4_arr, y_prev4_arr, h, f):
  y_n_1 = y_prev4_arr[-1] + h * (55/24 * f(t_prev4_arr[-1], y_prev4_arr[-1]) - 
                                      59/24 * f(t_prev4_arr[-2], y_prev4_arr[-2]) + 
                                      37/24 * f(t_prev4_arr[-3], y_prev4_arr[-3]) - 
                                      3/8 * f(t_prev4_arr[-4], y_prev4_arr[-4]))
  return y_n_1

def BrusselatorSolutionsAdams(t_prev4_arr, y_prev4_arr, N = int(1e4), h = 0.01):
  t_array = t_prev4_arr
  y_array = y_prev4_arr

  for i in range(3, N):
    y_new = AdamsIteration_4Order(t_prev4_arr, y_prev4_arr, h, Brusselator)
    y_array.append(y_new)
    t_array.append(h * (i + 1))

  return y_array, t_array

def ShowPlotsAdams(u, v, t):
  plt.figure(figsize=[12, 5], dpi=100)
  plt.plot(u, v)
  plt.title("Метод Адамса. Фазовый портрет v(u) при B = %.2f" %B_GLOBAL)
  plt.grid()
  plt.show()

  plt.figure(figsize=[12, 5], dpi=100)
  plt.plot(t, u, label = "u")
  plt.plot(t, v, label = "v")
  plt.xlabel("t")
  plt.ylabel("u, v")
  plt.legend()
  plt.title("Матод Адамса. u(t) и v(t) при B = %.2f" %B_GLOBAL)
  plt.grid()
  plt.show()

# ---------- main ----------

def main():
  N = 10000
  h = 0.01

  solution_arr, t_arr = BrusselatorSolutionsRK(N = N, h = h)
  solution_arr_transposed = np.array(solution_arr).transpose()
  ShowPlotsRK(solution_arr_transposed[0], solution_arr_transposed[1], t_arr)

  solution_arr, t_arr = BrusselatorSolutionsAdams(t_arr[0:4], solution_arr[0:4], N = N, h = h)
  solution_arr_transposed = np.array(solution_arr).transpose()
  ShowPlotsAdams(solution_arr_transposed[0], solution_arr_transposed[1], t_arr)


if __name__ == "__main__":
    main()