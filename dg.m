% Compute the global DG matrices Gs and Gf for the Legendre basis used on the problem
%   dt(U) + a*dx(U) + b*dx(P) = 0,
%   dt(P) + a*dx(P) + b*dx(U) = 0.
% The DG scheme is given in the form
%   v' = a/dx*Gs*v + b/dx*Gf*v
% where b>>a>0 and v=(u^T,p^T)^T are the degrees of freedom of (u,p).

% AUTHORS:
%   Slavko Brdar, University of Freiburg. 2013
%   Oswald Knoth, TROPOS, 2013

% LICENSE:
%   GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

clear 
clear global

global aadv = 0.6; % slow advection
global badv = 10.; % fast advection
global endt = 0.3; % model end time

global gridel = 128; % number of grid elements
global k = 2 % DG pol. degree

global tintegtype = 'imrk'
% (exrk) k=0:feuler, k=1:tvd2, k=2:tvd3, k=3:fe3, k=4:fe4, k=5:dp4
% (imrk) k=0:beuler, k=1:dirk22, k=2:dirk33
% (imexrk) k=1:imexspp2, k=2:imexspp3, k=3:imexmp, k=4:imexfb, k=5:imexd232
% (splitei) k=2:RK3
global tinteg = 'dirk33'
global kint = k;
global ns = 30;

global debug = 0; % standard check on initial data etc.
global deepdebug = 0; % even more debug (choose gridel small, i.e. 3)

global ndofs = gridel*(k+1); % total number of dofs for one variable
global dx = 1/gridel; % grid elements
global plotdx = 1/2000; % number of points used for plotting

n = (k+1)*gridel;

for i = 1:(k+1)
  for j = 1:(k+1)
    a(i,j) = 0.5 * (-1d0)^(j-1) * (2d0*i-1d0);
    bs(i,j) = -(2d0*i-1d0);
    bf(i,j) = 0.5 * (2d0*i-1d0)*(-1-(-1d0)^(j+i));
    cs(i,j) = (-1d0)^(i-1) * (2d0*i-1d0);
    e(i,j) = 0.5 * (2d0*i-1d0) * (-1d0+(-1d0)^(j+i));
    if ((i>j) && (mod(i-j,2)>0))
      e(i,j) = e(i,j) + (2d0*i-1d0)*2d0;
      bs(i,j) = bs(i,j) + (2d0*i-1d0)*2d0;
    end
  end
end

f = 0.5*cs;
cf = f;
d = -a;

lsmat=sparse(n,n);
mmat=sparse(n,n);
lfmat=sparse(n,n);
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
gsmat=sparse(2*n,2*n);
gfmat=sparse(2*n,2*n);
gsmat(1:n,1:n) = lsmat;

%lfmat=0;

%gsmat(n+1:2*n, 1:n) = 0;
%gsmat(1:n, n+1:2*n) = 0;
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


[u,p] = initialize(k,gridel,aadv,badv);

% plot numerical initial data after dg projection
plotdgsol(u,'dgu-init.png',k,gridel,plotdx);
figure(2)
%plotdgsol(p,'dgp-init.png',k,gridel,plotdx);


if (strcmp(tintegtype,'exrk'))
    cfl = 0.5;
    dt = cfl * dx / ((aadv+badv)*(2*k+1));
    [us,ps,dendt] = exrk(tinteg,kint,k,aadv,badv,dx,dt,endt,ndofs,u,p);
elseif (strcmp(tintegtype,'forback'))
    cfl = 0.5;
    dt = cfl * dx / ((aadv+badv)*(2*k+1));
    [us,ps,dendt] = forback(tinteg,kint,k,aadv,badv,dx,dt,endt,ndofs,u,p,gridel);
elseif (strcmp(tintegtype,'imrk'))
    cfl = 0.5;
    dt = 10. * cfl * dx / (max(aadv,1)*(2*k+1));
    [us,ps,dendt] = imrk(tinteg,kint,aadv,badv,dx,dt,endt,ndofs,u,p);
elseif (strcmp(tintegtype,'imexrk'))
    cfl = 0.5;
    dt = cfl * dx / (max(aadv,1)*(2*k+1));
    [us,ps,dendt] = imexrk(tinteg,kint,aadv,badv,dx,dt,endt,ndofs,u,p);
elseif (strcmp(tintegtype,'splitei'))
    cfl = 0.5;
    dt = cfl * dx / (max(aadv,1)*(2*k+1));
    [us,ps,dendt] = splitei(tinteg,kint,k,aadv,badv,dx,dt,endt,ndofs,u,p);
else
    printf('\nTime integration scheme \"%s\" is not available. Exiting...\n', tintegtype);
    exit(1);
end

% plot stability region for the time discretization scheme
%exrkplotstab(tinteg,k,dt,dx,aadv,badv);

% plot numerical initial data after dg projection
sum(abs(us))
plotdgsol(us,'dgu.png',k,gridel,plotdx);
hold on
plotdgsol(u,'dgu-init.png',k,gridel,plotdx);
%plotdgsol(ps,'dgp.png',k,gridel,plotdx);

% plot exact solution
fprintf('\nPlotting the analytical solution...');
x=0:plotdx:1;
exactsolt = @(x) exactsol(x,dendt,aadv,badv);
exactu = arrayfun(exactsolt,x);
plot(x,exactu);
exactsolname = 'exactu.png';
print(exactsolname,'-dpng');
fprintf('[ok] (see \"%s\")',exactsolname);

% compute error
err = l2err(us,k,dendt,aadv,badv,gridel);

fprintf('\n');
