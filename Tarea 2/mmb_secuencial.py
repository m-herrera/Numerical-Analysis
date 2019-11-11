import datetime
from numpy import *
from sympy import *
import numpy
from multiprocessing import *
import matplotlib.pyplot as plt


def subGradiente(gradiente, variables):
    solvers = []
    for i in range(0, len(variables)):
        sols = solve(gradiente[i], variables[i])
        solvers.append(sols)
    return solvers


def evalGradiente(gradiente, valores):
    """
    Esta funcion evalua en el vector gradiente los datos contenidos en el vector valores
     y retorna la norma de esta, esto con el objetivo de obtener el error de la funcion
    :param gradiente: Es un vector gradiente en forma simbolica
    :param variables: Son las variables simbolicas que forman parte de la gradiente
    :param valores: Son los valores en los cuales sera evaluado la gradiente
    :return: Es la norma del vector gradiente evaluada en valores
    """
    resultado = []
    for i in range(0, len(gradiente)):
        resultado.append(gradiente[i](*valores))
    resultado = numpy.array(resultado, dtype="float")
    return numpy.linalg.norm(resultado)


def calc_gradiente(funcion, variables):
    """
    Esta funcion es hecha para retornar el vector gradiente en forma simbolica
    :param funcion: LA funcion a calcular el gradiente
    :param variables: Son las variales al respecto a la que sera derivada
    :return: Un vector gradiente simbolico
    """
    symGradient = []
    gradient = []
    for i in range(0, len(variables)):
        partialDerivative = diff(funcion, variables[i])
        symGradient.append(partialDerivative)
        gradient.append(lambdify(variables, partialDerivative, "numpy"))
    return gradient, symGradient


def calcular_nuevoX(i, variables, f, vectorOG, solvers):
    """
    Este calcula un nuevo vector X pero solo optimizando la variable contenida en el espacio i
    mediante la regla de jacobi,tambien calcula el error de la nueva iteracion
    :param i: El numero de la variable en el vector que sera actualizada
    :param variable: Estas son las variables que seran reemplazadas con valores contenidos en vector
    :param funcion:Funcion a la cual se le esta calculando el minimo
    :param vector:El vector que contiene los datos  en los cuales sera evaluado la funcion
    :return: Una tupla con el error y la tupla
    """
    f = lambdify(variables, f, "numpy")
    vector = copy(vectorOG)
    soluciones = []
    for j in range(0, len(solvers[i])):
        funcion = lambdify(variables, solvers[i][j], "numpy")
        soluciones.append(funcion(*vector))
    valor = copy(vector)
    valor = numpy.array(valor, dtype='float')
    valor[i] = soluciones[0]
    for w in range(0, len(soluciones)):
        vector[i] = soluciones[w]
        if (f(*vector) < f(*valor)):
            valor[i] = soluciones[w]
    e = f(*valor)
    return (valor, e)

def mejoraMaxBloques(funcion, variables, vector, tol):
    """
    Este metodo calcula  el minimo de una funcion itilizando el metodo de la mejora maxima por bloques,
    pero esta funcion evalua cada una de los nuevos valores de forma secuencial, por lo cual es el metodo mas lento
    :param funcion: Funcion a la cual sra calculada el minimo
    :param variables: Un vector de variables simbolicas que forman parte de funcion
    :param vector: Un vector de valores iniciales que seran iterados
    :param tol: El valor en el cual se detiene la iteraciones
    :return:
    """

    g = calc_gradiente(funcion, variables)
    gradient = g[0]
    symGradient = g[1]
    solvers = subGradiente(symGradient, variables)
    ploterrors = []
    error = evalGradiente(gradient, vector)
    while ( error> tol):
        soluciones = []
        for i in range(0, len(vector)):
            soluciones.append(calcular_nuevoX(i, variables, funcion, vector, solvers))
        actual = soluciones[0]
        for i in range(0, len(soluciones)):
            if (soluciones[i][1] < actual[1]):
                actual = soluciones[i]
        ploterrors.append(error)
        vector = actual[0]
        error = evalGradiente(gradient, vector)

    print("Tiempo Finalizacion: ")
    print(datetime.datetime.now().time())
    plt.plot(ploterrors)
    plt.show()
    return vector

def ResolverNormal(N,tol):
    """
    Este es el metodo que resuelve el problema en cuestion de la tarea en forma secuencial
    por lo que crea la matriz y la crea
    :return:
    """
    variables = []
    for i in range(0, N):
        variables.append(Symbol("x" + str(i)))
    matriz = numpy.eye(N, dtype='object') * 6
    variablesMatriz = numpy.array(variables, dtype='object').T
    b = []
    for i in range(0, N):
        if (i == 0 or i == N-1):
            b.append(12)
        else:
            b.append(14)
    valoresparaeval = [1] * N
    for i in range(0, N-1):
        matriz[i][i + 1] = 2
    for i in range(0, N-1):
        matriz[i + 1][i] = 2
    bMatriz = numpy.array(b, dtype='object').T
    funcion = dot(dot(variablesMatriz.T, matriz), variablesMatriz)*1/2 - dot(bMatriz.T, variablesMatriz)
    print("Tiempo de inicio:")
    print(datetime.datetime.now().time())
    resultado = mejoraMaxBloques(funcion, variablesMatriz, valoresparaeval, tol)
    return resultado

"""
Aqui es donde se corre el codigo, se ejecuta la funcion ResolverNormal se escoje el tamaño de la matriz y tolerancia
"""
if __name__ == '__main__':
    print(ResolverNormal(50,0.01))
