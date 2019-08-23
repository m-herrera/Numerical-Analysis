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

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} funct3 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: kenne <kenne@DESKTOP-RK8LG59>
## Created: 2019-08-18

function [iteracion,aproximacion] = sne_fd_3 (funcion, tolerancia, valorInicial,graf=1)
  syms x;
    funcion = str2func( funcion);
    iteracion = 0;
    errores = [];
    iteraciones = [];
    while(abs(double(subs(funcion,x,valorInicial)))>=tolerancia);
       errores = [ errores double(abs(subs(funcion,x,valorInicial)))];
        iteraciones = [iteraciones iteracion];
        w= double(valorInicial + subs(funcion,x,valorInicial));
        divisor = double((subs(funcion,x,valorInicial) - subs(funcion,x,w))/(valorInicial-w));
        valorInicial =double( valorInicial - subs(funcion,x,valorInicial)/divisor);
        iteracion+=1;
    endwhile
    aproximacion = valorInicial;
    if(graf==1)
    plot(iteraciones,errores)
    xlabel('Number of iteration') 
    ylabel('Error') 
  endif
  if(!(graf==1 || graf ==0));
  print("WARNING: graf can only have two values, 0 or 1")
  endif
endfunction