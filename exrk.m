% explicit RK solver for system of ODEs
function [u,p,dendt] = exrk(tinteg,kInt,kPol,aadv,badv,dx,dt,endt,ndofs,un,pn)
  global gsmat gfmat plotdx gridel
  fprintf('\nExplicit \"%s\" Runge-Kutta time integration...',tinteg);
  gmat = aadv/dx*gsmat + badv/dx*gfmat;
  [a,b] = exbutcher(tinteg,kInt);
  nstages = max(1,size(b,2));

  checkexrk(tinteg,kInt);

  dendt = 0;
  tsteps = 0;
  vn = [un,pn];
  v = vn';
  iter=0;
  while (dendt-endt < -0.5*dt)
    iter=iter+1;
    kstage = zeros(size(vn,2),nstages);
    for stage = 1 : nstages
      for istage = 1 : stage-1
        kstage(:,stage) = kstage(:,stage) + a(stage,istage)*kstage(:,istage);
      end
      kstage(:,stage) = dt*kstage(:,stage) + v;
      kstage(:,stage) = gmat*kstage(:,stage);
    end
    v = v + dt*kstage*b';
    dendt = dendt + dt;
    tsteps = tsteps + 1;
    if mod(iter,50)==0
        dendt
        plotdgsol(v(1:ndofs)','dgu.png',kPol,gridel,plotdx);
    end
  end
  u = v(1:ndofs)';
  sum(abs(u))
  p = v(ndofs+1:end)';
  fprintf('[ok] (dt=%f, endt=%f, tsteps=%d)', dt, dendt, tsteps);
end



% Butcher tables
function [a,b] = exbutcher(tinteg,k)
  if (k == 0)
    if (strcmp(tinteg,'feuler'))
      a = [0];
      b = [1];
    end
  elseif (k == 1)
    if (strcmp(tinteg,'tvd2'))
      a = [0, 0; 1, 0];
      b = [1/2, 1/2];
    end
  elseif (k == 2)
    if (strcmp(tinteg,'tvd3'))
      a = [0, 0, 0; 1, 0, 0; 1/4, 1/4, 0];
      b = [1/6, 1/6, 2/3];
    end
  elseif (k == 3)
    if (strcmp(tinteg,'fe3'))
      % E. Fehlberg. Klassische Runge-Kutta Formeln vierter und niedrigerer Ordnung mit
      % Schrittweitenkontrolle und ihre Anwendung auf Warmeleitungsprobleme. Computing, 6:61-71, 1970
      a = [0, 0, 0, 0, 0, 0;...
           1/4, 0, 0, 0, 0, 0;...
           3/32, 9/32, 0, 0, 0, 0;...
           1932/2197, -7200/2197, 7296/2197, 0, 0, 0;...
           439/216, -8, 3680/513, -845/4104, 0, 0;...
           -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
      b = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    end
  elseif (k == 4)
    if (strcmp(tinteg,'fe4'))
      % E. Fehlberg. Klassische Runge-Kutta Formeln vierter und niedrigerer Ordnung mit
      % Schrittweitenkontrolle und ihre Anwendung auf Warmeleitungsprobleme. Computing, 6:61-71, 1970
      a = [0, 0, 0, 0, 0, 0;...
           1/4, 0, 0, 0, 0, 0;...
           3/32, 9/32, 0, 0, 0, 0;...
           1932/2197, -7200/2197, 7296/2197, 0, 0, 0;...
           439/216, -8, 3680/513, -845/4104, 0, 0;...
           -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
      b = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
    end
  elseif (k == 5)
    if (strcmp(tinteg,'dp4'))
      % J.R. Dormand, P.J. Prince. A family of embedded Runge-Kutta formulae. J. Comput. Appl. Math.,
      % 6:19-26, 1980
      a = [0, 0, 0, 0, 0, 0, 0;...
           1/5, 0, 0, 0, 0, 0, 0;...
           3/40, 9/40, 0, 0, 0, 0, 0;...
           44/45, -56/15, 32/9, 0, 0, 0, 0;...
           19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;...
           9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0;...
           35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
      b = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    end
  else
    printf('Explicit Runge-Kutta: butcher(%d) is not implemented', k);
    exit(1);
  end
end



function checkexrk(tinteg,k)
  % check if the Butcher table is strictly lower triangular
  [a,b] = exbutcher(tinteg,k);
  nstages = max(1,size(b,2));
  a
  for r = 1 : nstages
    for c = r : nstages
      if (abs(a(r,c))>1e-15)
        printf('exrk matrix is not strictly lower triangular! Exiting...');
        exit(1);
      end
    end
  end
end



% stability function for explicit Runge-Kutta
function [fr,fi] = exrkstab(tinteg,x,y,k,dt,dx,aadv,badv)
  [a,b] = exbutcher(tinteg,k);
  nstages = max(1,columns(b));

  checkexrk(tinteg,k);

  g  = (aadv+badv)/dx;
  kstage = zeros(1,nstages);
  for stage = 1 : nstages
    for istage = 1 : stage-1
      kstage(:,stage) = kstage(:,stage) + a(stage,istage)*kstage(:,istage);
    end
    kstage(:,stage) = (x+i*y)*dt*kstage(:,stage) + (x+y*i);
    kstage(:,stage) = g*kstage(:,stage);
  end
  fr = real(1 + dt*kstage*b');
  fi = imag(1 + dt*kstage*b');
end



% plot exrk stability
function exrkplotstab(tinteg,k,dt,dx,aadv,badv)
  printf('\nPlotting stability region for explicit Runge-Kutta...');
  scale = dx/(dt*(aadv+badv));
  [x,y] = meshgrid(-10*scale:0.075*scale:scale,-10*scale:0.075*scale:10*scale);
  f = @(x,y) exrkstab(tinteg,x,y,k,dt,dx,aadv,badv);
  [fr,fi] = arrayfun(f,x,y);
  g = fr.^2+fi.^2;
  titlename = strcat('Stability of explicit Runge-Kutta of order',' ',num2str(k));
  contour(x,y,g,[1 1],'b');
  title(titlename);
  xlabel('x');ylabel('y'); 
  filename = 'exrk-stab.png';
  print (filename,'-dpng');
  printf('[ok] (see \"%s\")', filename);
end
