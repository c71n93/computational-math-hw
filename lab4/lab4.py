import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

############### Imput Data ###############

years = np.array([1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000])
population = np.array([92228496, 106021537, 123202624, 132164569, 151325798, 179323175, 203211926, 226545805, 248709873, 281421906])

############### Newton interpolation ###############

def getPolinomCoefs(x, y):
    dividedDiffMatr = getPolinomDividedDiffs(x, y)
    size = dividedDiffMatr.shape[0]
    coefs = np.zeros(size)

    for i in range(size):
        coefs[i] = dividedDiffMatr[i][0]

    return coefs


def getPolinomDividedDiffs(x, y):
    dividedDiffMatr = np.zeros(shape = (y.size, y.size))

    for i in range(y.size): dividedDiffMatr[0][i] = y[i]

    for i in range(1, y.size):
        for j in range(y.size - i):
            dividedDiffMatr[i][j] = (dividedDiffMatr[i - 1][j + 1] - dividedDiffMatr[i - 1][j]) / (x[j + i] - x[j])
    return dividedDiffMatr


def polinom(coefs, x, point):
    result = 0

    for i in range(coefs.size):
        temp = coefs[i]
        for j in range(i):
            temp *= (point - x[j])
        result += temp

    return result


def drawNewtonsInterpolation(x, y):
    coefs = getPolinomCoefs(x, y)
    points = np.linspace(x[0], 2010, 1000)
    
    plt.figure(figsize=(6, 5))
    plt.title("Интерполяция Ньютона")
    plt.ylabel("$Численность$")
    plt.xlabel("$Год$")
    plt.plot(points, polinom(coefs, x, points), "black", label = "График полинома интерполяции")
    plt.plot(x, y, 'o', color='black', markersize = 4, label = "Исходные точки")
    plt.plot(2010, polinom(coefs, x, 2010), 'o', color='red', markersize = 4, label = "Экстраполированная точка")
    plt.legend(loc = "best")
    plt.grid()
    plt.savefig('images/newton.png')

    return


############### Spline interpolation ###############

def CalculateTriagCoeffsSpline(h, y):
    n = len(y)
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    f = np.zeros(n)

    a[0] = 0
    b[0] = 1
    c[0] = 0
    f[0] = 0

    for k in range(1, n - 1):
        a[k] = h[k] / 6
        b[k] = (h[k] + h[k + 1]) / 3
        c[k] = h[k + 1] / 6
        f[k] = (y[k + 1] - y[k]) / h[k + 1] - (y[k] - y[k - 1]) / h[k]

    a[n - 1] = 0
    b[n - 1] = 1
    c[n - 1] = 0
    f[n - 1] = 0

    return a, b, c, f


def SolveTriagSystem(a, b, c, f):
    n = len(a)
    p = np.zeros(n)
    r = np.zeros(n)

    p[0] = c[0] / b[0]
    r[0] = f[0] / b[0]

    for k in range(1, n - 1):
        p[k] = c[k] / (b[k] - a[k] * p[k - 1])
        r[k] = (f[k] - a[k] * r[k - 1]) / (b[k] - a[k] * p[k - 1])

    for k in range(n - 2, -1, -1):
        c[k] = r[k] - p[k] * c[k + 1]

    return c


def getSpline(x, y):
    n = len(x)
    h = np.zeros(n)

    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)

    for k in range(1, len(h)):
        h[k] = x[k] - x[k - 1]

    a_triag, b_triag, c_triag, f_triag = CalculateTriagCoeffsSpline(h, y)
    c = SolveTriagSystem(a_triag, b_triag, c_triag, f_triag)

    for k in range(1, n):
        a[k] = y[k]
        d[k] = (c[k] - c[k - 1]) / h[k]
        b[k] = 1/6 * (2 * c[k] + c[k - 1]) * h[k] + (y[k] - y[k - 1]) / h[k]


    def SplineFunc(s):
        k = 0
        if s <= x[1]:
            k = 1
        elif s > x[n - 1]:
            k = n - 1
        else:
            for i in range(1, n):
                if s > x[i - 1] and s <= x[i]:
                    k = i
                    break

        return a[k] + b[k] * (s - x[k]) + 0.5 * c[k] * ((s - x[k]) ** 2) + 1/6 * d[k] * ((s - x[k]) ** 3)

    return SplineFunc

def drawSpineInterpolation(x, y):
    spline = getSpline(x, y)
    points = np.linspace(x[0], 2010, 1000)

    plt.figure(figsize=(6, 5))
    plt.title("Сплайн-интерполяция")
    plt.ylabel("$Численность$")
    plt.xlabel("$Год$")
    plt.plot(points, [spline(points[i]) for i in range(points.size)], "black", label = "График полинома интерполяции")
    plt.plot(x, y, 'o', color='black', markersize = 4, label = "Исходные точки")
    plt.plot(2010, spline(2010), 'o', color='red', markersize = 4, label = "Экстраполированная точка")
    plt.legend(loc = "best")
    plt.grid()
    plt.savefig('images/spline.png')

    return

############### main ###############

def main():
    exactValue = 308745538

    coefs = getPolinomCoefs(years, population)
    extrapolatedNewton = polinom(coefs, years, 2010)

    spline = getSpline(years, population)
    extrapolatedSpline = spline(2010)

    drawNewtonsInterpolation(years, population)
    drawSpineInterpolation(years, population)

    print('Метод Ньютона: ', extrapolatedNewton)
    print('Сплайн-Интерполяция: ', extrapolatedSpline)
    print('Истинное значние: ', exactValue)

    return


if __name__ == "__main__":
    main()