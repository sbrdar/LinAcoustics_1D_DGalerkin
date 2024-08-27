% legendre polynomials on [-1,1]
function res = legendrefunc(x,n)
  if (n == 1)
    res = x;
  elseif (n == 0)
    res = 1;
  else
    res = x*(2-1/n).*legendrefunc(x,n-1) - (1-1/n)*legendrefunc(x,n-2);
  end
end
