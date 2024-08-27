% plot dg solution of order k on [0,1] for a given dof vector u
function plotdgsol(u,name,k,gridel,plotdx)
  fprintf('\nPlotting a dg solution...');
  x = 0:plotdx:1;
  dgsolval = @(x) evaldgsol(x,u,k,gridel);
  dgsol = arrayfun(dgsolval,x);
  plot(x,dgsol);
  set(gca,'xlim',[0,1],'ylim',[-1,1]);
  print(name,'-dpng');
  %fprintf('[ok] (see \"%s\")',name);
end

% evaluate dg solution at a given point
function dgval = evaldgsol(x,u,k,gridel)
  ge = min(gridel,max(1, ceil(x*gridel)));
  xl = (ge-1)/gridel;
  xr = ge/gridel;
  basis(1:k+1) = 0;
  for l = 1 : k+1
    basis(l) = legendrefunc(2/(xr-xl)*(x-(xl+xr)/2),l-1);
  end
  dgval = dot( u((ge-1)*(k+1)+1 : ge*(k+1)), basis );
end

% compute error in 'u' between our dg solution of order k and the exact analytical solution on [0,1]
function uerr = l2err(u,k,endt,aadv,badv,gridel)
  printf('\nComputing error...');
  uerr = 0;
  exactsolt = @(x) exactsol(x,endt,aadv,badv);
  dgsolt = @(x) evaldgsol(x,u,k,gridel);
  ul2err = @(x) (exactsolt(x)-dgsolt(x))^2;
  x = 0:0.01:1;
  y = arrayfun(ul2err,x);
  plot(x,y);
  print('err.png','-dpng');
  for ge = 1 : gridel
    xj = (ge-1)/gridel;
    xj1 = ge/gridel;
    [uerrint,uerrerr,nuerreval] = quad(ul2err,xj,xj1);
    uerr = uerr + uerrint;
  end
  uerr = sqrt(uerr);
  printf('[ok] (l2 error: %e)', uerr);
end
