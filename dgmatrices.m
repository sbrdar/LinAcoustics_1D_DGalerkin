% Compute the global DG matrices Gs and Gf for the Legendre basis used on the problem
%   dt(U) + a*dx(U) + b*dx(P) = 0,
%   dt(P) + a*dx(P) + b*dx(U) = 0.
% The DG scheme is given in the form
%   v' = a/dx*Gs*v + b/dx*Gf*v
% where b>>a>0 and v=(u^T,p^T)^T are the degrees of freedom of (u,p).

k
n = (k+1)*gridel;

for i = 1:(k+1)
  for j = 1:(k+1)
    a(i,j) = 0.5 * (-1d0)**(j-1) * (2d0*i-1d0);
    bs(i,j) = -(2d0*i-1d0);
    bf(i,j) = 0.5 * (2d0*i-1d0)*(-1-(-1d0)**(j+i));
    cs(i,j) = (-1d0)**(i-1) * (2d0*i-1d0);
    e(i,j) = 0.5 * (2d0*i-1d0) * (-1d0+(-1d0)**(j+i));
    if ((i>j) && (mod(i-j,2)>0))
      e(i,j) = e(i,j) + (2d0*i-1d0)*2d0;
      bs(i,j) = bs(i,j) + (2d0*i-1d0)*2d0;
    end
  end
end

f = 0.5*cs;
cf = f;
d = -a;

lsmat(1:n,1:n) = 0;
mmat(1:n,1:n) = 0;
lfmat(1:n,1:n) = 0;
for i = 1:gridel
  j = i;
  mmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = e;
  lfmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = bf;
  lsmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = bs;
  j = i+1;
  if (i == gridel)
    j = 1;
  end
  mmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = d;
  lfmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = a;
  j = i-1;
  if (i == 1)
    j = gridel;
  end
  mmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = f;
  lfmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = cf;
  lsmat((i-1)*(k+1)+1 : i*(k+1), (j-1)*(k+1)+1 : j*(k+1)) = cs;
end

% assemble global system matrices to be used later
global gsmat;
global gfmat;
gsmat(1:2*n,1:2*n) = 0;
gfmat(1:2*n,1:2*n) = 0;
gsmat(1:n,1:n) = lsmat;
gsmat(n+1:2*n, 1:n) = 0;
gsmat(1:n, n+1:2*n) = 0;
gsmat(n+1:2*n, n+1:2*n) = lsmat;
gfmat(1:n,1:n) = lfmat;
gfmat(n+1:2*n, 1:n) = mmat;
gfmat(1:n, n+1:2*n) = mmat;
gfmat(n+1:2*n, n+1:2*n) = lfmat;

if (deepdebug)
  printmatrix(a,'a');
  printmatrix(bs,'bs');
  printmatrix(bf,'bf');
  printmatrix(cf,'cf');
  printmatrix(cs,'cs');
  printmatrix(d,'d');
  printmatrix(e,'e');
  printmatrix(f,'f');
  printmatrix(mmat,'mmat');
  printmatrix(lfmat,'lfmat');
  printmatrix(lsmat,'lsmat');
  printmatrix(gsmat,'gsmat');
  printmatrix(gfmat,'gfmat');
  es = eig(lsmat);
  printmatrix(es,'es');
end

