% Butcher tables
function [a,b] = imbutcher(tinteg,k)
  if (k == 0)
    if (strcmp(tinteg,'beuler'))
      a = [1];
      b = [1];
    end
  elseif (k == 1)
    if (strcmp(tinteg,'dirk22'))
      alpha = 1+sqrt(2)/2;
      a = [alpha, 0; 1-alpha, alpha];
      b = [1-alpha, alpha];
    end
  elseif (k == 2)
    if (strcmp(tinteg,'dirk33'))
      % R. Alexander (Theorem 5). Diagonally implicit Runge-Kutta methods for stiff O.D.E.'s. SIAM J.
      % Numer. Anal., 14(6), pp. 1006-1021
      beta = 0.4358665214;
      tau = (1+beta)/2;
      c1 = -0.25*(6*beta*beta-16*beta+1);
      c2 = 0.25*(6*beta*beta-20*beta+5);
      a = [beta, 0, 0; tau-beta, beta, 0; c1, c2, beta];
      b = [c1, c2, beta];
    end
  else
    printf('Implicit Runge-Kutta: butcher(%d) is not implemented', k);
    exit(1);
  end
end




function checkimrk(tinteg,k)
  % check that we use DIRK
  [a,b] = imbutcher(tinteg,k);
  nstages = max(1,columns(b));
  for r = 1 : nstages
    if (a(1,1)!=a(r,r))
      printf('imrk matrix is not DIRK! Exiting...');
      exit(1);
    end
    for c = r+1 : nstages
      if (abs(a(r,c))>1e-15)
        printf('imrk matrix is not DIRK! Exiting...');
        exit(1);
      end
    end
  end
end




% stability function for implicit Runge-Kutta
function [fr,fi] = imrkstab(tinteg,x,y,k,dt,dx,aadv,badv)
  %global gsmat gfmat % too big to copy, access them directly -- do not change!
  %gmat = aadv/dx*gsmat + badv/dx*gfmat;
  [a,b] = imbutcher(tinteg,k);
  nstages = max(1,columns(b));

  g  = (aadv+badv)/dx;
  s = 1 - dt*a(1,1)*g;
  kstage = zeros(1,nstages);
  for stage = 1 : nstages
    for istage = 1 : stage-1
      kstage(:,stage) = (kstage(:,stage) + a(stage,istage)*kstage(:,istage)) / s;
    end
    kstage(:,stage) = (x+y*i)*dt*kstage(:,stage) + (x+y*i);
    kstage(:,stage) = g*kstage(:,stage);
  end
  fr = real(1 + dt*kstage*b');
  fi = imag(1 + dt*kstage*b');
end




% plot imrk stability
function imrkplotstab(tinteg,k,dt,dx,aadv,badv)
  printf('\nPlotting stability region for implicit Runge-Kutta...');

  checkimrk(tinteg,k);

  scale = dx/(dt*(aadv+badv));
  [x,y] = meshgrid(-10*scale:0.075*scale:scale,-10*scale:0.075*scale:10*scale);
  f = @(x,y) imrkstab(tinteg,x,y,k,dt,dx,aadv,badv);
  [fr,fi] = arrayfun(f,x,y);
  g = fr.^2+fi.^2;
  titlename = strcat('Stability of implicit Runge-Kutta of order',' ',num2str(k));
  contour(x,y,g,[1 1],'b');
  title(titlename);
  xlabel('x');ylabel('y'); 
  filename = 'imrk-stab.png';
  print ('imrk-stab.png','-dpng');
  printf('[ok] (see \"%s\")',filename);
end




% implicit RK solver for system of ODEs
function [u,p,dendt] = imrk(tinteg,kInt,aadv,badv,dx,dt,endt,ndofs,un,pn)
  global gsmat gfmat
  printf('\nImplicit \"%s\" Runge-Kutta time integration...',tinteg);
  gmat = aadv/dx*gsmat + badv/dx*gfmat;
  [a,b] = imbutcher(tinteg,kInt);
  nstages = max(1,columns(b));
  checkimrk(tinteg,kInt); % check that we use DIRK

  dendt = 0;
  tsteps = 0;
  vn = [un,pn];
  v = vn';
  s = eye(2*ndofs)-dt*a(1,1)*gmat;
  while (dendt - endt < -0.5*dt)
    kstage = zeros(columns(vn),nstages);
    for stage = 1 : nstages
      for istage = 1 : stage-1
        kstage(:,stage) = s \ (kstage(:,stage) + a(stage,istage)*kstage(:,istage));
      end
      kstage(:,stage) = dt*kstage(:,stage) + v;
      kstage(:,stage) = gmat*kstage(:,stage);
    end
    v = v + dt*kstage*b';
    dendt = dendt + dt;
    tsteps = tsteps + 1;
  end
  u = v'(1:ndofs);
  p = v'(ndofs+1:end);
  printf('[ok] (dt=%f, endt=%f, tsteps=%d)', dt, endt, tsteps);
end
