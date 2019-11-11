from sympy import *
from sympy import Symbol
import Functs


# print(wheneval)
#
# list = [1, 2, 3, 4, 5]
#
# i =
#
# newlist = list[:i] + list[i+1:]

# print(newlist)

variables = symbols("x y z")
eq = sympify("(x - 2) ** 4 + (x - 2 * y) ** 2 + z ** 2")
gradiente = Functs.calc_gradiente(eq, variables)[1]

# solutions = solve(gradiente[0], variables[0])
# pprint(solutions)
# wheneval = []
# for solution in solutions:
#     solution = lambdify(variables[1:], solution, "numpy")
#     wheneval.append(solution(25, 50))
# print(wheneval)


Functs.subGradiente(gradiente, variables, 1)

