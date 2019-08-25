from sympy import *
import matplotlib.pyplot as plt
import math

MAX_ITER = 1000


def set_max_iter(max_iter):
    global MAX_ITER
    MAX_ITER = max_iter


def sne_ud_1(x0, alpha, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        d2f = lambdify(variable, diff(funcion, variable, 2), "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            Lf = f(x) * d2f(x) / df(x)**2
            x = x - (1 + Lf / (2 * (1 - alpha * Lf))) * f(x) / df(x)
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
            plot(errors, "Halley's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


# Chun's Method
def sne_ud_2(x0, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            y = x - f(x) / df(x)
            x = y - (f(x) + 2 * f(y)) / f(x) * f(y) / df(x)
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
            plot(errors, "Chun's Method")
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_ud_3(x0, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            y = x - f(x) / df(x)
            x = y - (f(x) + f(y))**2 / (f(x)**2 - 5 * f(y)**2) * f(y) / df(x)
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
            plot(errors, "Title")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_ud_4(x0, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            if df(x) == 0:
                print("Error, the derivative undefined the function, f'(x)==0\n")
                break
            y = x - f(x) / df(x)
            x = y - f(y) / df(x)
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
            plot(errors, "Title")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_ud_5(x0, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]
        while error >= tol and iteracion < MAX_ITER:
            denominador = df(x - (1 / 2) * f(x) / df(x))
            x = x - f(x) / denominador
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
            plot(errors, "Title")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion


#Weerakoon and Fernando method
def sne_ud_6(x0, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]
        while error > tol and iteracion < MAX_ITER:
            denominador = (df(x) + df(x - f(x)/df(x)))
            x = x - (2 * f(x)) / denominador
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
            plot(errors, "Weerakoon and Fernando")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion



# Dong's method
# Requiere conocer la multiplicidad (m) de la raiz
def sne_ud_7(x0, m, tol, funcion, graf = 1):
    iteracion = 0
    x = x0
    try:
        if tol < 0:
            print("Error: Tolerance value must be positive\n")
            return
        variable = Symbol("x")
        funcion = sympify(funcion)
        f = lambdify(variable, funcion, "numpy")
        df = lambdify(variable, diff(funcion, variable), "numpy")
        error = abs(f(x))
        errors = [error]

        if m <= 0:
            print("ERROR: m must be greater than 0")
            return None

        while error > tol and iteracion < MAX_ITER:
            denominador = df(x)
            if denominador == 0:
                break
            y = x - math.sqrt(m) * f(x) / denominador
            x = y - m * ((1 - 1 / math.sqrt(m)) ** (1 - m)) * (f(y) / denominador)
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
            plot(errors, "Dong")
        except:
            print("Unable to plot errors")
    return x, iteracion


def plot(errors, title):
    plt.plot(errors)
    plt.title(title)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()

