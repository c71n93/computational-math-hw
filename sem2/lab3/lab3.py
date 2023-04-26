import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
import time

def System(x, vec):
  f = np.array([0.0, 0.0])
  f[0] = vec[0]
  f[1] = -0.5 * vec[1]**2 / (1 - 0.5 * vec[0])
  return f

def RungeKuttaIteration_4Order(x_n, vec_n, f, h = 0.01):
  f_1 = f(x_n, vec_n)
  f_2 = f(x_n + h / 2, vec_n + f_1 * h / 2)
  f_3 = f(x_n + h / 2, vec_n + f_2 * h / 2)
  f_4 = f(x_n + h, vec_n + f_3 * h)

  vec_n_1 = vec_n + h / 6 * (f_1 + 2 * f_2 + 2 * f_3 + f_4)
  return vec_n_1

def SolutionsRK(vec_start, x_start = 0.0, x_end = 1.0, h = 0.01):
  x_num = int(x_end / h)

  x_array = np.array([x_start])
  vec_array = np.array([vec_start])

  for i in range(0, x_num):
    vec_new = RungeKuttaIteration_4Order(x_array[i], vec_array[i], System, h)
    vec_array = np.append(vec_array, [vec_new], axis = 0)
    x_array = np.append(x_array, [h * (i + 1)])

  return vec_array, x_array

def BinSolutionSearch(z0_left, z0_right, y0, y1 = 0, abs_tol = 0.001):
  z0_middle = (z0_left + z0_right) / 2
  new_sol, new_x_array = SolutionsRK(vec_start = np.array([y0, z0_middle]))
  y1_new = new_sol[-1][0] # y(1)
  print("y'(0) = " + str(z0_middle) + " y(1) = " + str(y1_new))
  if abs(y1_new - y1) < abs_tol:
    return new_sol, new_x_array, z0_middle
  elif y1_new > y1:
    return BinSolutionSearch(z0_middle, z0_right, y0)
  else:
    return BinSolutionSearch(z0_left, z0_middle, y0)
  
def IterationSolutionSearch():
  y0 = 0.25
  arr = np.arange(-1, 0.05, 0.05)
  for i in arr:
    new_sol, new_x_array = SolutionsRK(vec_start = np.array([y0, i]))
    y_array = new_sol.transpose()[0]
    ShowSoulutionPlots(y_array, new_x_array, new_sol[0][1])
    print("")
    
def ShowSoulutionPlots(y, x, p):
  # clear_output(wait=True)
  plt.figure(figsize=[12, 5], dpi=100)
  plt.plot(x, y, label = "y")
  plt.xlabel("x")
  plt.ylabel("y")
  plt.legend()
  plt.title("Метод стрельбы. y(x), y'(0) = " + str(p))
  plt.grid()
  plt.show()

# y0 = 0.25
# sol_array, x_array, p = BinSolutionSearch(-1.0, 0.0, y0 = y0)
# y_array = sol_array.transpose()[0]
# ShowSoulutionPlots(y_array, x_array, p)
IterationSolutionSearch()