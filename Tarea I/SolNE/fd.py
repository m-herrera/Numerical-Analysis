from sympy import *
from Plotter import *
from sympy.parsing.sympy_parser import parse_expr


# improved Ostrowski’s method free from derivatives
def sne_fd_1(x0, tol, function, graf):
    variable = Symbol("x")
    function = sympify(function)
    f = lambdify(variable, function, "numpy")
    iteration = 0
    x = x0
    error = abs(f(x))
    errors = [error]
    while error >= tol:
        y = x - (2 * f(x) ** 2) / (f(x + f(x)) - f(x - f(x)))
        z = y - ((y - x) * f(y)) / (2 * f(y) - f(x))
        x = z - ((y - x) * f(z)) / (2 * f(y) - f(x))
        iteration += 1
        error = abs(f(x))
        errors.append(error)
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0")
    elif graf:
        plot(errors, "Improved Ostrowski's Method")
    return x, iteration


# steffensen's method
def sne_fd_2(x0, tol, function, graf):
    variable = Symbol("x")
    function = sympify(function)
    f = lambdify(variable, function, "numpy")
    iteration = 0
    x = x0
    error = abs(f(x))
    errors = [error]
    while error >= tol:
        x = x - (f(x) ** 2) / (f(x + f(x)) - f(x))
        iteration += 1
        error = abs(f(x))
        errors.append(error)
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0")
    elif graf:
        plot(errors, "Steffensen's Method")
    return x, iteration


def sne_fd_3(valorInicial, tolerancia,funcion,graf = 1):
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


def sne_fd_4(valorInicial,anterior, tolerancia,funcion,graf = 1):
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
def sne_fd_5(x0, tol, funcion, graf=1):
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
def sne_fd_6(x0, gamma, tol, funcion, graf=1):
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
