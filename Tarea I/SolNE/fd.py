import math
import sympy
from sympy import *
import matplotlib.pyplot as plt
from sympy.parsing.sympy_parser import parse_expr
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


def sne_fd_4(valorInicial, tolerancia,funcion,graf = 1):
    x = Symbol('x')
    funcion = parse_expr( funcion)
    iteracion = 0
    puntos = []
    errors = []
    while(abs(N(funcion.subs(x,valorInicial)))>=tolerancia):
        w= N(valorInicial + funcion.subs(x,valorInicial))
        divisor = N((funcion.subs(x,valorInicial) - funcion.subs(x,w))/(valorInicial-w))
        valorInicial =N( valorInicial - funcion.subs(x,valorInicial)/divisor)
        iteracion+=1
        errors.append(abs(funcion.subs(x,valorInicial)))
        puntos.append(valorInicial)
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(errors,puntos)
        

def sne_fd_5(valorInicial,anterior, tolerancia,funcion,graf = 1):
    x = Symbol('x')
    funcion = parse_expr( funcion)
    iteracion = 0
    puntos = []
    errors = []
    while(abs(N(funcion.subs(x,valorInicial)))>=tolerancia):
        divisor = (funcion.subs(x,anterior)- funcion.subs(x,2*valorInicial-anterior))/(anterior - (2*valorInicial-anterior))
        calculo =N( valorInicial - funcion.subs(x,valorInicial)/divisor)
        anterior = valorInicial
        valorInicial = calculo
        iteracion+=1
        errors.append(abs(funcion.subs(x,valorInicial)))
        puntos.append(valorInicial)
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(errors,puntos)


# Jain's method
def sne_fd_6(x0, tol, funcion, graf=1):
    variable = Symbol("x")
    funcion = sympify(funcion)
    iteracion = 0
    f = lambdify(variable, funcion, "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]
    while error > tol and iteracion < 1000:
        denominadory = (f(x + f(x)) - f(x))
        if denominadory == 0:
            break;
        y = x - (f(x) ** 2) / denominadory
        denominadorx = ((f(x + f(x)) - f(x)) * (f(x) - f(y)))
        if denominadorx == 0:
            break;
        x = x - (f(x) ** 3) / denominadorx
        iteracion += 1
        error = abs(f(x))
        errors.append(error)
    if (graf != 0 and graf != 1):
        print("WARNING: graf has two possible values, 1 or 0")
    elif (graf == 1):
        plot(errors, "Jain")
    return x, iteracion

# Zhang's method
def sne_fd_7(x0, gamma, tol, funcion, graf=1):
    variable = Symbol("x")
    funcion = sympify(funcion)
    iteracion = 0
    f = lambdify(variable, funcion, "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]

    def divDif1(x, y):
        return (f(x) - f(y)) / (x - y)

    def divDif2(x, y, z):
        return (divDif1(x, y) - divDif1(y, z)) / (x - z)

    while error > tol and iteracion < 1000:
        w = x + gamma * f(x)
        denominadory = divDif1(x, w)
        if denominadory == 0:
            break;
        y = x - f(x) / denominadory
        denominadorx = divDif1(x, y) + (y - x) * divDif2(x, w, y)
        if denominadorx == 0:
            break;
        x = y - f(y) / denominadorx

        iteracion += 1
        error = abs(f(x))
        errors.append(error)
    if (graf != 0 and graf != 1):
        print("WARNING: graf has two possible values, 1 or 0")
    elif (graf == 1):
        plot(errors, "Zhang")
    return x, iteracion



# print(sne_fd_3(-0.5, 0.000001, "lambda x: x**2 - math.e**x-3*x+2"))
#print(sne_fd_3(2, 0.00000000000001, "lambda x: math.atan(x)"))
def plot_f(errors,values):
    plt.plot(values, errors)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()
    
def plot(errors, title):
    plt.plot(errors)
    plt.title(title)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()