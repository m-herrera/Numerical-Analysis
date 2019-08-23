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

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} sne_ud_1 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: mherr <mherr@MSI>
## Created: 2019-08-22

# Halley's Method
function [x, iteration] = sne_ud_1(x0, alpha, tol, f, graf)
    syms x;
    f = str2func(f);
    df = matlabFunction(diff(f, x));
    d2f = matlabFunction(diff(f, x, 2));
    x = x0;
    iteration = 0;
    error = abs(f(x));
    errors = [error];
    while error >= tol
        Lf = f(x) * d2f(x)/ df(x)^2;
        x = x - (1 + Lf / (2 * (1 - alpha * Lf))) * f(x) / df(x);
        iteration += 1;
        error = abs(f(x));
        errors = horzcat(errors, error);
    endwhile
    if graf != 0 && graf != 1
        disp("WARNING: graf has two possible values, 1 or 0");
    elseif graf
      figure
      plot(errors);
      title("UD1")
      xlabel('Number of iteration') 
      ylabel('Error') 
    end
endfunction