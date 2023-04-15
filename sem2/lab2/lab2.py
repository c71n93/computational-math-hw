import numpy as np
from matplotlib import pyplot as plt

# Task conditions
epsilon = 0.0001
x_0     = 1.0
y_0     = 5.0
alpha_0 = 0.001
T_k     = 1500

# Rosenbrok method constants
p_1  = 0.435866521508459
p_2  = 0.4782408332745185
p_3  = 0.0858926452170225
a    = p_1
b_21 = p_1
b_31 = p_1
b_32 = -2.116053335949811


def System(vec):
  x, y, alpha = vec
  dxdt = x * (1 - 0.5 * x - (2 / (7 * alpha**2)) * y)
  dydt = y * (2 * alpha - 3.5 * alpha**2 * x - 0.5 * y)
  dadt = epsilon * (2 - 7 * alpha * x)
  return np.array([dxdt, dydt, dadt])


def Jac(vec):
    x, y, alpha = vec
    J = np.array([
        [1 - x - (2 / (7 * alpha**2)) * y, (-2 / (7 * alpha**2)) * x, (-4 / (7 * alpha**3)) * x * y],
        [-3.5 * alpha**2 * y, 2 * alpha - 3.5 * alpha**2 * x - 0.5 * y, -0.5 * y],
        [-7 * epsilon * x, epsilon * (2 - 7 * alpha), -7 * epsilon * alpha * x]
    ])
    return J


def Dn(vec, h):
    return np.eye(3) + a * h * Jac(vec)


def Rosenbrok3Order(start_vec, h, n_of_steps):
    res = np.zeros((n_of_steps, np.shape(start_vec)[0]))
    res[0] = start_vec
    
    for n in range(0, n_of_steps - 1):
        D_n = Dn(res[n], h)
        k_1 = np.linalg.solve(D_n, h * System(res[n]))
        k_2 = np.linalg.solve(D_n, h * System(res[n] + b_21 * k_1))
        k_3 = np.linalg.solve(D_n, h * System(res[n] + b_31 * k_1 + b_32 * k_2))
        res[n + 1] = res[n] + p_1 * k_1 + p_2 * k_2 + p_3 * k_3
    return res


def showPlotRosenbrok(res, t):
  res = res.T
  plt.figure(figsize=[16, 10])
  plt.xlabel("t")
  plt.ylabel("x")
  plt.minorticks_on()
  plt.tight_layout()
  labels = ["x(t)", "y(t)", "a(t)"]
  for i in range(len(res)):
    plt.plot(t, res[i], label=labels[i])
  plt.legend(loc = 'best', fontsize = 12)
  plt.show()


def showSinglePlotRosenbrok(single_res, t, label):
  plt.figure(figsize=[16, 10])
  plt.xlabel("t")
  plt.ylabel(label)
  plt.minorticks_on()
  plt.tight_layout()
  plt.plot(t, single_res, label=label)
  plt.legend(loc = 'best', fontsize = 12)
  plt.show()


def main():
  h = 0.05
  n_of_steps = int(T_k / h) + 1
  t = np.linspace(0, T_k, n_of_steps)

  start_vec = np.array([x_0, y_0, alpha_0])
  res = Rosenbrok3Order(start_vec, h, n_of_steps)
  showPlotRosenbrok(res, t)
  # showSinglePlotRosenbrok(res.T[0], t, "a")

if __name__ == "__main__":
    main()
