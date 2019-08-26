from SolNE.ud import *
from SolNE.fd import *

print(sne_ud_1(100, 0, 0.0001, "E**x+3*x-sin(x)", 1))
print(sne_ud_2(100, 0.0001, "E**x+3*x-sin(x)", 1))
print(sne_ud_3(100, 0.0001, "E**x+3*x-sin(x)", 1))
print(sne_ud_4(100, 0.0001, "E**x+3*x-sin(x)", 1))
print(sne_ud_5(100, 0.0001, "E**x+3*x-sin(x)", 1))
print(sne_ud_6(100, 1, 0.0001, "E**x+3*x-sin(x)", 1))

print(sne_ud_1(5, 1/2, 0.0001, "x**x + tan(x)", 1))
print(sne_ud_2(5, 0.0001, "x**x + tan(x)", 1))
print(sne_ud_3(5, 0.0001, "x**x + tan(x)", 1))
print(sne_ud_4(5, 0.0001, "x**x + tan(x)", 1))
print(sne_ud_5(5, 0.0001, "x**x + tan(x)", 1))
print(sne_ud_6(5, 1, 0.0001, "x**x + tan(x)", 1))

print(sne_ud_1(25, 0, 0.0001, "ln(x**2) + x - 5", 1))
print(sne_ud_2(25, 0.0001, "ln(x**2) + x - 5", 1))
print(sne_ud_3(25, 0.0001, "ln(x**2) + x - 5", 1))
print(sne_ud_4(25, 0.0001, "ln(x**2) + x - 5", 1))
print(sne_ud_5(25, 0.0001, "ln(x**2) + x - 5", 1))
print(sne_ud_6(25, 1, 0.0001, "ln(x**2) + x - 5", 1))

print(sne_fd_1(100, 0.000000001, "x**2+3*x-10", 1))
print(sne_fd_2(100, 0.000000001, "x**2+3*x-10", 1))
#print(sne_fd_3(100, 0.000000001, "x**2+3*x-10", 1))
print(sne_fd_4(100, 0, 0.000000001, "x**2+3*x-10", 1))
print(sne_fd_5(100, 0.000000001, "x**2+3*x-10", 1))
print(sne_fd_6(100, 1, 0.000000001, "x**2+3*x-10", 1))

print(sne_fd_1(1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
print(sne_fd_2(1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
#print(sne_fd_3(1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
print(sne_fd_4(2, 1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
print(sne_fd_5(1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
print(sne_fd_6(1, 1, 0.000000001, "x**4+sin(pi/x**2)-5", 1))
