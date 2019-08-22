import sympy
from sympy import *
import matplotlib.pyplot as plt
from sympy.parsing.sympy_parser import parse_expr
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

def sne_ud_4(valorInicial, tolerancia,funcion,graf = 1):
    x = Symbol('x')
    funcion = parse_expr( funcion)
    iteracion = 0
    puntos = []
    errors = []
    derivada = diff(funcion, x)
    while(abs(funcion.subs(x,valorInicial))>=tolerancia):
        if(derivada.subs(x,valorInicial) == 0):
            print("Error, the derivative undifined the function, f'(x)==0"
                  )
        y = valorInicial - funcion.subs(x,valorInicial)/derivada.subs(x,valorInicial)
        valorInicial = N(y - funcion.subs(x,y)/derivada.subs(x,valorInicial))
        puntos.append( valorInicial)
        errors.append(abs(funcion.subs(x,valorInicial)))
        iteracion+=1
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(errors,puntos)


def sne_ud_5(valorInicial, tolerancia,funcion,graf = 1):
    x = Symbol('x')
    funcion = parse_expr( funcion)
    iteracion = 0
    puntos = []
    errors = []
    derivada = diff(funcion, x)
    while(abs(N(funcion.subs(x,valorInicial)))>=tolerancia):
        if(derivada.subs(x,valorInicial) == 0):
            print("Error, the derivative undifined the function, f'(x)==0")
        divisor = N(derivada.subs(x,valorInicial-(1/2) * funcion.subs(x,valorInicial)/derivada.subs(x,valorInicial)))
        valorInicial = N(valorInicial - funcion.subs(x,valorInicial)/divisor)
        puntos.append(valorInicial)
        errors.append(abs(funcion.subs(x,valorInicial)))
        iteracion+=1
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(errors,puntos)
        
def plot_f(errors,values):
    plt.plot(values, errors)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()
