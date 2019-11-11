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
## @deftypefn {} {@var{retval} =} thomas (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: kenne <kenne@DESKTOP-RK8LG59>
## Created: 2019-11-10

function solucion = thomas (matriz, vector,largo)
p = 0;
q = 0;
P = [];
Q = [];
for i = 1 : largo
  if(i > 1)
    a = matriz(i, i - 1);
  else
    a = 0;
  endif
  b = matriz(i, i);
  
  q =  (vector(i) - q * a) / (b - p * a);
  Q = [Q, q];
  
  if i != rows(matriz)
    c = matriz(i, i + 1);
    p = c / (b - p * a);
    P = [P, p];
  else
    P = [P, 0];
  endif    
endfor
 
x = [];
xi = 0;
for i = largo:-1:1
  xi = Q(i) - P(i) * xi;
  x = [xi, x];
endfor
solucion=x

endfunction
