function [x, iteracion] = sne_fd_5(x0, tol, f, graf=1)
  f = str2func(f);
  pkg load symbolic;
  syms x;
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];
  while error > tol && iteracion < 1000
    denominadory = f(x + f(x)) - f(x);
    if denominadory == 0
      break;
    endif;
    y = x - (f(x) ** 2) / denominadory;
    denominadorx = (f(x + f(x)) - f(x)) * (f(x) - f(y));
    if denominadorx == 0
      break;
    endif;
    x = x - (f(x) ** 3) / denominadorx;
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
  title("Jain method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction