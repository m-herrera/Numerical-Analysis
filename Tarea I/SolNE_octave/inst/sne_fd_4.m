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

##Implementación del método mejorado de Kurchatov
## Entradas:
## -x0: valor inicial de iteración (tipo: numérico real)
## -prev: valor inicial de iteración (tipo: numérico real)
## -tol: tolerancia mínima del error (tipo: numérico real positivo)
## -funcion: función sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
## -graf: bandera para graficar o no el error de la función (tipo: numérico 1 o 0)
## Salidas:
## -Valor aproximado de la solución según tolerancia indicada o hasta que se indefina el procedimiento
## -Cantidad de iteraciones realizadas según tolerancia indicada o hasta que se indefina el procedimiento
## -Gráfica del error en función del número de iteraciones
function [x, iteracion] = sne_fd_4(x0, prev, tol, funcion, graf=1)
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
  
    if tol < 0
      disp("Error: Tolerance value must be positive")
      return
    endif
  while error > tol && iteracion < 1000
    denominador = (f(prev) - f(2 * x - prev))/(prev - (2 * x - prev));
    temp = x - f(x) / denominador;
    prev = x;
    x = temp;
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
  title("Kurchatov method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction
