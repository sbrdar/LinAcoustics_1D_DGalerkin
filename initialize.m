% initialize the global vector of dofs v=(u^T,p^T)^T
function [u,p] = initialize(k,gridel,aadv,badv)
  fprintf('\nProjecting intial data...');
  for ge = 1 : gridel
    xj = (ge-1)/gridel;
    xj1 = xj + 1/gridel;
    for l = 1 : (k+1)
      ufunc = @(x) uinit(x,badv).*legendrefunc(2/(xj1-xj).*(x-(xj+xj1)/2),l-1);
      pfunc = @(x) pinit(x,badv).*legendrefunc(2/(xj1-xj).*(x-(xj+xj1)/2),l-1);
      uint = quad( ufunc, xj, xj1 );
      pint = quad( pfunc, xj, xj1 );
      u((ge-1)*(k+1)+l) = (2*l-1)/2*uint*2*gridel;
      p((ge-1)*(k+1)+l) = (2*l-1)/2*pint*2*gridel;
    end
  end
  fprintf('[ok] (pol.deg.=%d)',k);
end

% u initial data
function u = uinit(x,badv)
  k = 14*pi/0.1;
  x1 = 0.25;
  x0 = 0.75;
  s0 = 0.1;
  p0 = exp(-((x-x0)./s0).^2);
  p1 = exp(-((x-x1)./s0).^2) .* cos(x*k);
  p = p0 + p1;
  u = p;
end

% p initial data
function p = pinit(x,badv)
  k = 14*pi/0.1;
  x1 = 0.25;
  x0 = 0.75;
  s0 = 0.1;
  p0 = exp(-((x-x0)./s0).^2);
  p1 = exp(-((x-x1)./s0).^2) .* cos(x*k);
  p = p0 + p1;
end

function [u,p] = exactsol(x,t,aadv,badv)
  lambda1 = aadv - badv;
  lambda2 = aadv + badv;
  x1 = x-t*lambda1-floor(x-t*lambda1);
  x2 = x-t*lambda2-floor(x-t*lambda2);
  uow = @(x) -uinit(x,badv)+pinit(x,badv);
  pow = @(x) uinit(x,badv)+pinit(x,badv);
  u = 0.5*(-uow(x1)+pow(x2));
  p = 0.5*(uow(x1)+pow(x2));
end



% if (debug)
%   printf('\nPlotting the analytical initial data...');
%   x=0:plotdx:1;
%   uinitb = @(x) uinit(x,badv)
%   yu = arrayfun(uinitb,x);
%   plot(x,yu);
%   print('u-init.png','-dpng');
%   pinitb = @(x) pinit(x,badv)
%   yp = arrayfun('pinitb',x);
%   plot(x,yp);
%   print('p-init.png','-dpng');
%   printf('[ok]');
% end
