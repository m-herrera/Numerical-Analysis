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

## Implementaci�n del m�todo de Traub
## Entradas:
## -x0: valor inicial de iteraci�n (tipo: num�rico real)
## -tol: tolerancia m�nima del error (tipo: num�rico real positivo)
## -funcion: funci�n sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
## -graf: bandera para graficar o no el error de la funci�n (tipo: num�rico 1 o 0)
## Salidas:
## -Valor aproximado de la soluci�n seg�n tolerancia indicada o hasta que se indefina el procedimiento
## -Cantidad de iteraciones realizadas seg�n tolerancia indicada o hasta que se indefina el procedimiento
## -Gr�fica del error en funci�n del n�mero de iteraciones

function [x, iteracion] = sne_ud_3(x0, tol, f, graf=1)
  try
    f = str2func(f);
  catch
    disp("Error: syntax error in the function handle.")
    disp(lasterr)
    x = 0;
    iteration = 0;
    return
  end_try_catch
  try
  pkg load symbolic;
  syms x;
  derivative = diff(f, x);
  df = matlabFunction(derivative);
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];

  if tol < 0
    disp("Error: Tolerance value must be positive")
    return
  endif
  while error > tol && iteracion < 1000
    if df(x) == 0
        print("Error, the derivative undefined the function, f'(x)==0\n")
        break
    endif
    y = x - f(x) / df(x);
    x = y - f(y) / df(x);
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
  title("Traub method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction
