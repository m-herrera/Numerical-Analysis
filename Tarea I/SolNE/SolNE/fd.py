from sympy import *
import matplotlib.pyplot as plt

MAX_ITER = 1000


def set_max_iter(max_iter):
    global MAX_ITER
    MAX_ITER = max_iter


# improved Ostrowski’s method free from derivatives
def sne_fd_1(x0, tol, funcion, graf=1):
    """Implementación del método mejorado de Ostrowski libre de derivadas
    Entradas:
    -x0: valor inicial de iteración (tipo: numérico real)
    -tol: tolerancia mínima del error (tipo: numérico real positivo)
    -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
    -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
    Salidas:
    -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
    -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
    -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            y = x - (2 * f(x) ** 2) / (f(x + f(x)) - f(x - f(x)))
            z = y - ((y - x) * f(y)) / (2 * f(y) - f(x))
            x = z - ((y - x) * f(z)) / (2 * f(y) - f(x))
            iteracion += 1
            error = abs(f(x))
            errors.append(error)
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Improved Ostrowski's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


# steffensen's method
def sne_fd_2(x0, tol, funcion, graf=1):
    """Implementación del método de Steffensen libre de derivadas
        Entradas:
        -x0: valor inicial de iteración (tipo: numérico real)
        -tol: tolerancia mínima del error (tipo: numérico real positivo)
        -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
        -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
        Salidas:
        -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
        -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
        -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            x = x - (f(x) ** 2) / (f(x + f(x)) - f(x))
            iteracion += 1
            error = abs(f(x))
            errors.append(error)
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Steffensen's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_fd_3(x0, a, tol, funcion, graf=1):
    """Implementación del método mejorado de Ren
            Entradas:
            -x0: valor inicial de iteración (tipo: numérico real)
            -a: parámetro numérico
            -tol: tolerancia mínima del error (tipo: numérico real positivo)
            -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
            -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
            Salidas:
            -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
            -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
            -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]

        def divDif1(x, y):
            return (f(x) - f(y)) / (x - y)

        while error >= tol and iteracion < MAX_ITER:
            y = x - (f(x) ** 2) / (f(x + f(x)) - f(x))
            z = x - (2 * f(x) ** 2) / (f(x + f(x)) - f(x - f(x)))
            x = y - f(y) / (divDif1(x, y) + divDif1(y, z) - divDif1(x, z) + a * (y - x) * (y - z))
            iteracion += 1
            error = abs(f(x))
            errors.append(error)
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Ren's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_fd_4(x0, prev, tol, funcion, graf = 1):
    """Implementación del método mejorado de Kurchatov
        Entradas:
        -x0: valor inicial de iteración (tipo: numérico real)
        -prev: valor inicial de iteración (tipo: numérico real)
        -tol: tolerancia mínima del error (tipo: numérico real positivo)
        -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
        -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
        Salidas:
        -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
        -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
        -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            denominador = (f(prev) - f(2 * x - prev))/(prev - (2 * x - prev))
            temp = x - f(x) / denominador
            prev = x
            x = temp
            error = abs(f(x))
            errors.append(error)
            iteracion += 1
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Kurchatov's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


# Jain's method
def sne_fd_5(x0, tol, funcion, graf = 1):
    """Implementación del método de Jain
        Entradas:
        -x0: valor inicial de iteración (tipo: numérico real)
        -tol: tolerancia mínima del error (tipo: numérico real positivo)
        -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
        -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
        Salidas:
        -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
        -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
        -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]
        while error > tol and iteracion < MAX_ITER:
            denominadory = f(x + f(x)) - f(x)
            y = x - (f(x) ** 2) / denominadory
            denominadorx = (f(x + f(x)) - f(x)) * (f(x) - f(y))
            x = x - (f(x) ** 3) / denominadorx
            iteracion += 1
            error = abs(f(x))
            errors.append(error)
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Jain's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


# Zhang's method
def sne_fd_6(x0, gamma, tol, funcion, graf = 1):
    """Implementación del método de Zhang
        Entradas:
        -x0: valor inicial de iteración (tipo: numérico real)
        -gamma: parametro numérico
        -tol: tolerancia mínima del error (tipo: numérico real positivo)
        -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
        -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
        Salidas:
        -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
        -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
        -Gráfica del error en función del número de iteraciones"""
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        error = abs(f(x))
        errors = [error]

        def divDif1(x, y):
            return (f(x) - f(y)) / (x - y)

        def divDif2(x, y, z):
            return (divDif1(x, y) - divDif1(y, z)) / (x - z)

        while error > tol and iteracion < MAX_ITER:
            w = x + gamma * f(x)
            denominadory = divDif1(x, w)
            y = x - f(x) / denominadory
            denominadorx = divDif1(x, y) + (y - x) * divDif2(x, w, y)
            x = y - f(y) / denominadorx
            iteracion += 1
            error = abs(f(x))
            errors.append(error)
    except ZeroDivisionError:
        print("Error: Division by zero\nShowing partial result\n")
    except OverflowError:
        print("Error: Iteration overflows due to initial value being too large\nShowing partial result\n")
    except:
        print("Error: invalid input\nShowing partial result\n")
    if graf != 0 and graf != 1:
        print("WARNING: graf has two possible values, 1 or 0\n")
    elif graf:
        try:
            plot(errors, "Zhang's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


def plot(errors, title):
    plt.plot(errors)
    plt.title(title)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()

