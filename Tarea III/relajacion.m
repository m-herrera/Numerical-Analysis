function x = relajacion(A, b, w = -1, x0 = -1, tol)
  n = length(A);
  
  # valor por defecto (0, 0, 0, 0)
  if (x0 == -1)
    x0 = zeros(n,1);
  endif
  
  # Obtener descomposición L, D, U
  L = D = U = zeros(n);
  for i = 1:n
    for j = 1:n
      if i > j
        L(i, j) = A(i, j);
      elseif i < j
        U(i, j) = A(i, j);
      else
        D(i, j) = A(i, j);
      endif
    endfor
  endfor
  
  
  if (w == -1)
    # Obtener autovalores (Av)
    [~, Av] = eig(inv(D)*(-(L+U)));
    # Calcular w
    w = 2/(1 + sqrt(1 - max(Av(:))**2))
  endif
  
  M_inv = inv(L + D/w);
  N = ((1 - w) / w * D - U);
  
  x = x0;
  error = norm(A * x - b);
  
  # iterar 
  while error >= tol
    x = M_inv * N * x + M_inv * b;
    error = norm(A * x - b);
  endwhile
  
endfunction
