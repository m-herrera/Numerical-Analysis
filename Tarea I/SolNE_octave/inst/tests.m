## Copyright (C) 2019 jassonrm
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.
##
## Author: jassonrm <MBP-de-Jasson>
## Created: 2019-08-18

## Prueba del paquete computacional resolviendo una simple ecuacion cuadratica

[x1, k1] = sne_ud_1(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)
[x2, k2] = sne_ud_2(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x3, k3] = sne_ud_3(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x4, k4] = sne_ud_4(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x5, k5] = sne_ud_5(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x6, k6] = sne_ud_6(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)


[x1, k1] = sne_fd_1(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x2, k2] = sne_fd_2(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x3, k3] = sne_fd_3(100, 2, 0.0001, "@(x)x**2+3*x-10", 1)
[x4, k4] = sne_fd_4(100, 101, 0.0001, "@(x)x**2+3*x-10", 1)
[x5, k5] = sne_fd_5(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x6, k6] = sne_fd_6(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)