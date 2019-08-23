function [x, iteracion] = sne_ud_6(x0, tol, f, graf=1)
  f = str2func(f);
  pkg load symbolic;
  syms x;
  derivative = diff(f, x);
  df = matlabFunction(derivative);
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];
  while error > tol && iteracion < 1000
    denominadorx = (df(x) + df(x - f(x)/df(x)));
    if denominadorx == 0
      break;
    endif;
    x = x - (2 * f(x))/ denominadorx;
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
  title("Weerakoon and Fernando method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction