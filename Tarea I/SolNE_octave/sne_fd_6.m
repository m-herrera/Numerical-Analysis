function [x, iteracion] = sne_fd_6(x0, gamma, tol, f, graf=1)
  f = str2func(f);
  pkg load symbolic;
  syms x;
  iteracion = 0;
  x = x0;
  error = abs(f(x));
  errors = [error];
  
  function [result] = divDif1(a, b)
    result = (f(a) - f(b)) / (a - b);
  endfunction

  function [result] = divDif2(a, b, c)
    result = (divDif1(a, b) - divDif1(b, c)) / (a - c);
  endfunction

  
  while error > tol && iteracion < 1000
    w = x + gamma * f(x);
    denominadory = divDif1(x, w);
    if denominadory == 0
      break;
    endif;
    y = x - f(x) / denominadory;
    denominadorx = divDif1(x, y) + (y - x) * divDif2(x, w, y);
    if denominadorx == 0
      break;
    endif;
    x = y - f(y) / denominadorx;
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
  title("Zhang method")
  xlabel('Number of iteration') 
  ylabel('Error') 
endfunction
