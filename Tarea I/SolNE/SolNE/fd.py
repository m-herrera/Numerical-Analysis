from sympy import *
import matplotlib.pyplot as plt

MAX_ITER = 1000


def set_max_iter(max_iter):
    global MAX_ITER
    MAX_ITER = max_iter


# improved Ostrowski’s method free from derivatives
def sne_fd_1(x0, tol, funcion, graf=1):
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


def sne_fd_3(x0, tol, funcion, graf=1):
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
            w = x + f(x)
            denominador = (f(x) - f(w)) / (x - w)
            x = x - f(x) / denominador
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
            plot(errors, "title")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion


def sne_fd_4(x0, prev, tol, funcion, graf = 1):
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
            temp = x - f(x0) / denominador
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
            plot(errors, "title")  # TODO
        except:
            print("Unable to plot errors")
    return x, iteracion


# Jain's method
def sne_fd_5(x0, tol, funcion, graf = 1):
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
            plot(errors, "Jain")
        except:
            print("Unable to plot errors")
    return x, iteracion


# Zhang's method
def sne_fd_6(x0, gamma, tol, funcion, graf = 1):
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
            plot(errors, "Zhang")
        except:
            print("Unable to plot errors")
    return x, iteracion


def plot(errors, title):
    plt.plot(errors)
    plt.title(title)
    plt.ylabel("Error")
    plt.xlabel("Iteraciones")
    plt.show()
