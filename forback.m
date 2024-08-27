% explicit RK solver for system of ODEs
function [u,p,dendt] = forback(tinteg,kPol,aadv,badv,dx,dt,endt,ndofs,un,pn,gridel)
  global gsmat gfmat plotdx % too big to copy, access them directly -- do not change!
  fprintf('\nExplicit \"%s\" Runge-Kutta time integration...',tinteg);
  gmat = aadv/dx*gsmat + badv/dx*gfmat;
  
  dendt = 0;
  tsteps = 0;
  uu = un';
  pp = pn';
  iter = 0;
  while (dendt-endt < -0.5*dt)
      iter=iter+1;
      uu = uu + dt*gmat(1:ndofs,:)*[uu;pp];
      pp = pp + dt*gmat(ndofs+1:2*ndofs,:)*[uu;pp];
      
    dendt = dendt + dt;
    tsteps = tsteps + 1;
    if mod(iter,500)==0
        plotdgsol(uu','dgu.png',kPol,gridel,plotdx);
    end
  end
  u = uu';
  p = pp'';
  fprintf('[ok] (dt=%f, endt=%f, tsteps=%d)', dt, dendt, tsteps);
end

