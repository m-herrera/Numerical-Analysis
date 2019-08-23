function [x, iteracion] = sne_ud_7(x0, m, tol, f, graf=1)
  f = str2func(f);
  pkg load symbolic;
  syms x;
  derivative = diff(f, x);
  df = matlabFunction(derivative);
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];
  
  if m <= 0
    print("ERROR: m must be greater than 0")
  endif
  
  while error > tol && iteracion < 1000
    denominador = df(x);
    if denominador == 0
      break;
    endif;
    y = x - sqrt(m) * f(x) / denominador;
    x = y - m * ((1 - 1 / sqrt(m)) ** (1 - m)) * (f(y) / denominador);
    iteracion += 1;
    error = abs(f(x));
    errors = [errors, error];
  endwhile
  if graf != 1 && graf != 0
    print("WARNING: graf has two possible values, 1 or 0")
  elseif graf == 0
    return
  endif
  figure
  plot(errors);
  title("Dong method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction