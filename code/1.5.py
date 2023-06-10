import numpy
from scipy import integrate


def f(x):
    return 1.0 / (2.0 + x)


def rectangular_method(a, b):
    n = b - a  # количество отрезков (при h = 1)
    h = 1
    I = 0  # интегральная сумма
    x = a + h / 2
    for i in range(n):
        x += h  # левая граница
        I += f(x)  # ((b - a) / n) * f(x(i-1) + h / 2)

    return I

def trapezoidal_method(a, b):
    h = 1
    n = b - a // h  # количество отрезков (при h = 1)
    I = 0
    for i in range(n):
        x1 = a + h * i
        x2 = x1 + h
        I += (f(x1) + f(x2)) / 2 * (x2 - x1)
    return I

def find_para(xi1, xi2, xi3):
    a = 0
    b = 0
    c = 0
    x = numpy.array([xi1, xi2, xi3])
    y = numpy.array([f(xi1), f(xi2), f(xi3)])
    z = numpy.polyfit(x, y, 2)
    return z[0], z[1], z[2]


def integrte_para(a, b, c, start, end):
    return ((a * end ** 3) / 3 + (b * end ** 2) / 2 + c * end) - \
        ((a * start ** 3) / 3 + (b * start ** 2) / 2 + c * start)


def simpson_method(a, b):
    I = 0
    h = 1
    n = (b - a) // (2 * h)
    for i in range(n):
        xi1 = a + 2 * h * i
        xi2 = xi1 + h
        xi3 = xi1 + 2 * h
        a_p, b_p, c_p = find_para(xi1, xi2, xi3)
        I += integrte_para(a_p, b_p, c_p, xi1, xi3)
    return I


def weddle(f, a, b, n):
    h = (b - a) / n
    s = 0.0
    x1 = a
    x2 = a + h
    while (x2 <= b):
        dx = (x2 - x1) / 6
        s = s + (h / 20) * (f(x1) + 5 * f(x1 + dx) + f(x1 + 2 * dx) + \
        6 * f(x1 + 3 * dx) + f(x1 + 4 * dx) + 5 * f(x1 + 5 * dx) + f(x1 + 6 * dx))
        x1 = x2
        x2 = x1 + h
    return s

def find_para_7(xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8):
    x = numpy.array([xi1, xi2, xi3, xi4, xi5, xi5, xi6, xi7, xi8])
    y = numpy.array([f(xi1), f(xi2), f(xi3), f(xi4), f(xi5), f(xi5), f(xi6), f(xi7), f(xi8)])
    z = numpy.polyfit(x, y, 7)
    return z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]


def integrte_para_7(a1, a2, a3, a4, a5, a6, a7, a8, start, end):
    return ((a1 * end ** 8) / 8 + (a2 * end ** 7) / 7 + (a3 * end ** 6) / 6 + (a4 * end ** 5) / 5 + (
                a5 * end ** 4) / 4 + (a6 * end ** 3) / 3 + (a7 * end ** 2) / 2 + (a8 * end ** 1) / 1) - (
            (a1 * start ** 8) / 8 + (a2 * start ** 7) / 7 + (a3 * start ** 6) / 6 + (a4 * start ** 5) / 5 + (
            a5 * start ** 4) / 4 + (a6 * start ** 3) / 3 + (a7 * start ** 2) / 2 + (a8 * start ** 1) / 1
    )


def simpson_method_7(a, b):
    I = 0
    h = 1
    # n = (b - a) // (2 * h)
    n = 2
    h = (b - a) / 14
    for i in range(n):
        xi1 = a + 7 * h * i
        xi2 = xi1 + h
        xi3 = xi1 + 2 * h
        xi4 = xi1 + 3 * h
        xi5 = xi1 + 4 * h
        xi6 = xi1 + 5 * h
        xi7 = xi1 + 6 * h
        xi8 = xi1 + 7 * h
        a1, a2, a3, a4, a5, a6, a7, a8 = find_para_7(xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8)
        I += integrte_para_7(a1, a2, a3, a4, a5, a6, a7, a8, xi1, xi7)
    return I

a = -1
b = 3
result = rectangular_method(a, b)
predicted = integrate.quad(f, a, b)[0]

print("{} rectangular_method || fault {}".format(result, predicted - result))

result = trapezoidal_method(a, b)

print("{} trapezoidal_method || fault {}".format(result, predicted - result))

result = integrate.quad(f, a, b)

print("{} library".format(result))

result = simpson_method(a, b)

print("{} simpson_method || fault {}".format(result, predicted - result))

result = weddle(f, a, b, 4)

print("{} weddle || fault {}".format(result, predicted - result))
