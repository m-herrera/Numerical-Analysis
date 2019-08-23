from sympy import *
from Plotter import *
from sympy.parsing.sympy_parser import parse_expr


# Halley's Method
def sne_ud_1(x0, alpha, tol, function, graf):
    variable = Symbol("x")
    function = sympify(function)
    iteration = 0
    f = lambdify(variable, function, "numpy")
    df = lambdify(variable, diff(function, variable), "numpy")
    d2f = lambdify(variable, diff(function, variable, 2), "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]
    while error >= tol:
        Lf = f(x) * d2f(x) / df(x)**2
        x = x - (1 + Lf / (2 * (1 - alpha * Lf))) * f(x) / df(x)
        iteration += 1
        error = abs(f(x))
        errors.append(error)
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0")
    elif graf:
        plot(errors, "Halley's Method")
    return x, iteration


# Chun's Method
def sne_ud_2(x0, tol, function, graf):
    variable = Symbol("x")
    function = sympify(function)
    iteration = 0
    f = lambdify(variable, function, "numpy")
    df = lambdify(variable, diff(function, variable), "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]
    while error >= tol:
        y = x - f(x) / df(x)
        x = y - (f(x) + 2 * f(y)) / f(x) * f(y) / df(x)
        iteration += 1
        error = abs(f(x))
        errors.append(error)
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0")
    elif graf:
        plot(errors, "Chun's Method")
    return x, iteration


def sne_ud_3(x0, tol, function, graf):
    variable = Symbol("x")
    function = sympify(function)
    iteration = 0
    f = lambdify(variable, function, "numpy")
    df = lambdify(variable, diff(function, variable), "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]
    while error >= tol:
        y = x - f(x) / df(x)
        x = y - (f(x) + f(y))**2 / (f(x)**2 - 5 * f(y)**2) * f(y) / df(x)
        iteration += 1
        error = abs(f(x))
        errors.append(error)
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0")
    elif graf:
        plot(errors, "Nameless Method")
    return x, iteration

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
        iteraciones.append(abs(funcion.subs(x,valorInicial)))
        iteracion+=1
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(iteraciones,puntos)


def sne_ud_5(valorInicial, tolerancia,funcion,graf = 1):
    x = Symbol('x')
    funcion = parse_expr( funcion)
    iteracion = 0
    puntos = []
    iteraciones = []
    derivada = diff(funcion, x)
    while(abs(N(funcion.subs(x,valorInicial)))>=tolerancia):
        if(derivada.subs(x,valorInicial) == 0):
            print("Error, the derivative undifined the function, f'(x)==0")
        divisor = N(derivada.subs(x,valorInicial-(1/2) * funcion.subs(x,valorInicial)/derivada.subs(x,valorInicial)))
        valorInicial = N(valorInicial - funcion.subs(x,valorInicial)/divisor)
        puntos.append(valorInicial)
        iteraciones.append(iteracion)
        iteracion+=1
    print("The number of iterations is: " + str(iteracion))
    print("With an aproximation of: " +str(valorInicial))
    if(graf!= 0 and graf!=1):
        print("WARNING: graf has two possible values, 1 or 0")
    if(graf==1):
        plot_f(iteraciones, puntos)


#Funciones Jasson
#Weerakoon and Fernando method
def sne_ud_6(x0, tol, funcion, graf=1):
    variable = Symbol("x")
    funcion = sympify(funcion)
    iteracion = 0
    f = lambdify(variable, funcion, "numpy")
    df = lambdify(variable, diff(funcion, variable), "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]
    while error > tol and iteracion < 1000:
        denominadorx = (df(x) + df(x - f(x)/df(x)))
        if denominadorx == 0:
            break
        x = x - (2 * f(x))/ denominadorx
        iteracion += 1
        error = abs(f(x))
        errors.append(error)

    if (graf != 0 and graf != 1):
        print("WARNING: graf has two possible values, 1 or 0")
    elif (graf == 1):
        plot(errors, "Weerakoon and Fernando")
    return x, iteracion


# Dong's method
#Requiere conocer la multiplicidad (m) de la raiz
def sne_ud_7(x0, m, tol, funcion, graf=1):
    variable = Symbol("x")
    funcion = sympify(funcion)
    iteracion = 0
    f = lambdify(variable, funcion, "numpy")
    df = lambdify(variable, diff(funcion, variable), "numpy")
    x = x0
    error = abs(f(x))
    errors = [error]

    if m <= 0:
        print("ERROR: m must be greater than 0")
        return None

    while error > tol and iteracion < 1000:
        denominador = df(x)
        if denominador == 0:
            break
        y = x - math.sqrt(m) * f(x) / denominador
        x = y - m * ((1 - 1 / math.sqrt(m)) ** (1 - m))  * (f(y) / denominador)
        iteracion += 1
        error = abs(f(x))
        errors.append(error)
    if (graf != 0 and graf != 1):
        print("WARNING: graf has two possible values, 1 or 0")
    elif (graf == 1):
        plot(errors, "Dong")
    return x, iteracion
