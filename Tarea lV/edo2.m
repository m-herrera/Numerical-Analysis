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
## @deftypefn {} {@var{retval} =} edo2 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: kenne <kenne@DESKTOP-RK8LG59>
## Created: 2019-11-10
#edo2("@(x) -1/x","@(x) 1/(4*x**2)-1","@(x) 0",0.1,1,6,1,6)
function [x,y] = edo2 (p,q,f,h,a,b,y0,yn)
  p = str2func(p);
  q = str2func(q);
  f = str2func(f);
  vectorx= [a];
  respaldoA= a;
  while(a+h<b)
    a= a+h;
    vectorx = [vectorx a];
  endwhile
  largo=length (vectorx);
  e0 = (h/2*p(vectorx(1)) +1)*y0;
  en = (-1*h/2*p(vectorx(largo)) +1)*yn;
  matriz = zeros(largo);
  x=vectorx;
  for i=(1:largo)
    matriz(i,i)= 2+h**2*q(vectorx(i));
  endfor
  for i = (2:largo)
    matriz(i,i-1)= (h/-2) * p(vectorx(i)) -1;
  endfor
  for i = (1:largo-1)
    matriz(i,i+1) =  (h/2)*p(vectorx(i)) -1;
  endfor
  matriz=matriz;
  y=[];
  resultado = -1*h**2*f(vectorx(1)) +e0;
  y=[y resultado];
  for i=2:largo-2
    y= [ y -1*h**2*f(vectorx(i))];
  endfor
  y= [y -1*h**2*f(vectorx(largo)) +en ];
  y= thomas(matriz,y,largo);
endfunction
