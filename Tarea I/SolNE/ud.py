import sympy

#Halley's Method
def sne_ud_1(x0, tol, func):
    x = sympy.symbols('x')
    x_k = x0
    deriv1 = sympy.diff(eval(func)(x), x)
    deriv2 = sympy.diff(deriv1, x)
    n_iter = 0
    a=0
    while abs(eval(func)(x_k)) >= tol:
        Lf = eval(func)(x_k)*deriv2.evalf(subs={x: x_k})/(deriv1.evalf(subs={x: x_k})**2)
        x_k = x_k - (1 + 0.5 * Lf / (1 - a*Lf)) * eval(func)(x_k) / deriv1.evalf(subs={x: x_k})
        n_iter += 1
    return x_k, n_iter

#Chun's Method
def sne_ud_2(x0, tol, func):
    x = sympy.symbols('x')
    x_k = x0
    deriv1 = sympy.diff(eval(func)(x), x)
    n_iter = 0
    while abs(eval(func)(x_k)) >= tol:
        y_k = x_k - eval(func)(x_k) / deriv1.evalf(subs={x: x_k})
        x_k = y_k - (eval(func)(x_k) + 2*eval(func)(y_k))/eval(func)(x_k)*eval(func)(y_k) / deriv1.evalf(subs={x: x_k})
        n_iter += 1
    return x_k, n_iter


def sne_ud_3(x0, tol, func):
    x = sympy.symbols('x')
    x_k = x0
    deriv1 = sympy.diff(eval(func)(x), x)
    n_iter = 0
    while abs(eval(func)(x_k)) >= tol:
        y_k = x_k - eval(func)(x_k) / deriv1.evalf(subs={x: x_k})
        x_k = y_k - (eval(func)(x_k) + eval(func)(y_k))**2/(eval(func)(x_k)**2-5*eval(func)(y_k)**2) * eval(func)(y_k) / deriv1.evalf(subs={x: x_k})
        n_iter += 1
    return x_k, n_iter

