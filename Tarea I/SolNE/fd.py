import math


# improved Ostrowski’s method free from derivatives
def sne_fd_1(x0, tol, func):
    x_k = x0
    n_iter = 0
    while abs(eval(func)(x_k)) >= tol:
        y_k = x_k - (2 * eval(func)(x_k) ** 2) / (eval(func)(x_k + eval(func)(x_k)) - eval(func)(x_k - eval(func)(x_k)))
        z_k = y_k - ((y_k - x_k) * eval(func)(y_k)) / (2 * eval(func)(y_k) - eval(func)(x_k))
        x_k = z_k - ((y_k - x_k) * eval(func)(z_k)) / (2 * eval(func)(y_k) - eval(func)(x_k))
        n_iter += 1
    return x_k, n_iter

#Ostrowski’s method free from derivatives
def sne_fd_2(x0, tol, func):
    x_k = x0
    n_iter = 0
    while abs(eval(func)(x_k)) >= tol:
        y_k = x_k - (2 * eval(func)(x_k) ** 2) / (eval(func)(x_k + eval(func)(x_k)) - eval(func)(x_k - eval(func)(x_k)))
        x_k = x_k - (2 * eval(func)(x_k) ** 2) / (eval(func)(x_k + eval(func)(x_k)) - eval(func)(x_k - eval(func)(x_k))) * (eval(func)(y_k)-eval(func)(x_k))/(2*eval(func)(y_k)-eval(func)(x_k))
        n_iter += 1
    return x_k, n_iter

#steffensen's method
def sne_fd_3(x0, tol, func):
    x_k = x0
    n_iter = 0
    while abs(eval(func)(x_k)) >= tol:
        x_k = x_k - (eval(func)(x_k) ** 2) / (eval(func)(x_k + eval(func)(x_k)) - eval(func)(x_k))
        n_iter += 1
    return x_k, n_iter



print(sne_fd_3(-0.5, 0.000001, "lambda x: x**2 - math.e**x-3*x+2"))
#print(sne_fd_3(2, 0.00000000000001, "lambda x: math.atan(x)"))
