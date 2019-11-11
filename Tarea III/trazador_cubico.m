function S = trazador_cubico(x, y)
  n = length(x);
  ti = x;
  x = y;
  h = ti(2) - ti(1);
  
  A = zeros(n);
  for i = 1: n
    for j = 1: n
      if i == 1 || i == n
        if i == j
          A(i, j) = 1;
        endif
      elseif i == j
            A(i, j) = 4;
      elseif i == j + 1 || j == i + 1
            A(i, j) = 1;
      endif
    endfor
  endfor
  
  k = 4;
  b = zeros(n, 1);
  for i = 1: n
    if i == 1 || i == n
      b(i) = 0;
    else
      b(i) = 6/h**2 *(x(k-3) - 2*x(k-2) + x(k-1));
      k += 1;
    endif
  endfor
  
  S = 0;
  z = relajacion(A, b, :, :, 0.001);
  t = sym('t');
  for i = 1: n - 1
    S = S + (t>vpa(ti(i)))* (t<=vpa(ti(i+1)))*(vpa(z(i+1)/(6*h))*(t-vpa(ti(i)))**3 + vpa(z(i)/(6*h))*(vpa(ti(i+1))-t)**3 + vpa((x(i+1)/h - h*z(i+1)/6))*(t-vpa(ti(i))) + vpa((x(i)/h - h*z(i)/6))*(vpa(ti(i+1)) - t));
  endfor
  
endfunction
