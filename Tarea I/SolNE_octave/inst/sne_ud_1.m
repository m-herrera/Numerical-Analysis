## Copyright (C) 2019 mherr
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
## Author: mherr <mherr@MSI>
## Created: 2019-08-22

## Implementaci�n del m�todo de Halley y variaciones
## Entradas:
## -x0: valor inicial de iteraci�n (tipo: num�rico real)
## -alpha: selecciona entre variaciones del metodo, 0 para Chebyshev's, 1/2 para Halley's y 1 para Super-Halley's
## -tol: tolerancia m�nima del error (tipo: num�rico real positivo)
## -funcion: funci�n sobre la cual iterar, siguiendo lineamientos de sympy (tipo: cadena de caracteres)
## -graf: bandera para graficar o no el error de la funci�n (tipo: num�rico 1 o 0)
## Salidas:
## -Valor aproximado de la soluci�n seg�n tolerancia indicada o hasta que se indefina el procedimiento
## -Cantidad de iteraciones realizadas seg�n tolerancia indicada o hasta que se indefina el procedimiento
## -Gr�fica del error en funci�n del n�mero de iteraciones
function [x, iteration] = sne_ud_1(x0, alpha, tol, funcion, graf)
    pkg load symbolic;
    syms x;
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
    pkg load symbolic;
    df = matlabFunction(diff(f, x));
    d2f = matlabFunction(diff(f, x, 2));
    x = x0;
    iteration = 0;
    error = abs(f(x));
    errors = [error];
    if tol < 0
      disp("Error: Tolerance value must be positive")
      return
    endif
    while error >= tol
        Lf = f(x) * d2f(x)/ df(x)^2;
        x = x - (1 + Lf / (2 * (1 - alpha * Lf))) * f(x) / df(x);
        iteration += 1;
        error = abs(f(x));
        errors = horzcat(errors, error);
    endwhile
    catch
      disp("WARNING: math error. Showing partial result")
      disp(lasterr)
    end_try_catch
    if graf != 0 && graf != 1
        disp("WARNING: graf has two possible values, 1 or 0");
    elseif graf
      figure
      plot(errors);
      title("Halley method")
      xlabel('Number of iteration') 
      ylabel('Error') 
    end
endfunction