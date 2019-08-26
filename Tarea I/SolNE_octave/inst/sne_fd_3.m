## Copyright (C) 2019 kenne
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
## Author: kenne <kenne@DESKTOP-RK8LG59>
## Created: 2019-08-18

## Implementación del método mejorado de Ren
## Entradas:
## -x0: valor inicial de iteración (tipo: numérico real)
## -a: parametro
## -tol: tolerancia mínima del error (tipo: numérico real positivo)
## -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
## -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
## Salidas:
## -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
## -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
## -Gráfica del error en función del número de iteraciones
function [x, iteracion] = sne_fd_3(x0, a, tol, funcion, graf=1)
  try
    f = str2func(funcion);
  catch
    disp("Error: syntax error in the function handle.")
    disp(lasterr)
    x = 0;
    iteration = 0;
    return
  end_try_catch
  try
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];
  
  function [result] = divDif1(a, b)
    result = (f(a) - f(b)) / (a - b);
  endfunction
  
    if tol < 0
      disp("Error: Tolerance value must be positive")
      return
    endif
  while error > tol && iteracion < 1000
    y = x - (f(x) ** 2) / (f(x + f(x)) - f(x));
    z = x - (2 * f(x) ** 2) / (f(x + f(x)) - f(x - f(x)));
    x = y - f(y) / (divDif1(x, y) + divDif1(y, z) - divDif1(x, z) + a * (y - x) * (y - z));
    iteracion += 1;
    error = abs(f(x));
    errors = [errors, error];
  endwhile
      catch
      disp("WARNING: math error. Showing partial result")
      disp(lasterr)
    end_try_catch
  if graf != 1 && graf != 0
    print("WARNING: graf has two possible values, 1 or 0")
  elseif graf == 0
    return
  endif
  figure
  plot(errors);
  title("Ren method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction