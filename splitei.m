% explicit RK solver for system of ODEs
function [u,p,dendt] = splitei(tinteg,kInt,kPol,aadv,badv,dx,dt,endt,ndofs,un,pn)
global gsmat gfmat plotdx gridel ns % too big to copy, access them directly -- do not change!
fprintf('\nSplitExplicit \"%s\" Runge-Kutta time integration...',tinteg);
%gmat = aadv/dx*gsmat + badv/dx*gfmat;
[alpha,beta,gamma,d] = spbutcher(tinteg,kPol-1);
nstages = max(1,size(beta,2))-1;

%checkexrk(tinteg,kInt);

dendt = 0;
tsteps = 0;
yn = [un,pn];
Y=zeros(2*ndofs,nstages+1);
FY=zeros(2*ndofs,nstages);
F=zeros(2*ndofs);
iter=0;
while (dendt-endt < -0.5*dt)
    for stage=1:nstages+1
        Y(:,stage)=yn';
    end
    for stage=1:nstages
        FY(:,stage)=aadv/dx*gsmat*Y(:,stage);
        for istage=2:stage
            Y(:,stage+1)=Y(:,stage+1)+alpha(stage+1,istage)*(Y(:,istage)-yn');
        end
        F=0.0e0;
        for istage=2:stage
            F=F+gamma(stage+1,istage)/dt*(Y(:,istage)-yn');
        end
        for istage=1:stage
            F=F+beta(stage+1,istage)*FY(:,istage);
        end
        if d(stage+1)~=0.0e0
            F=F/d(istage+1);
            nsFB=ceil(ns*d(stage+1));
            dTau=dt*d(stage+1)/nsFB;
            [Y(1:ndofs,stage+1), Y(ndofs+1:2*ndofs,stage+1)] =...
                forback1(Y(1:ndofs,stage+1)',Y(ndofs+1:2*ndofs,stage+1)'...
               ,F(1:ndofs)',F(ndofs+1:2*ndofs)'...
               ,badv,dx,dTau,nsFB);
        end
    end 
    yn=Y(:,nstages+1)';
    dendt = dendt + dt;
    tsteps = tsteps + 1;
    if mod(tsteps,50)==0
        dendt
        plotdgsol(yn(1:ndofs)','dgu.png',kPol,gridel,plotdx);
    end

end
dendt
plotdgsol(yn(1:ndofs)','dgu.png',kPol,gridel,plotdx);
u = Y(1:ndofs,nstages)';
p = Y(ndofs+1:end,nstages)';
end

% Butcher tables
function [alpha,beta,gamma,d] = spbutcher(tinteg,k)
switch tinteg
    case{'RK3'}
       d=[   0.0000000000000000       0.33333333333333331       0.50000000000000000        1.0000000000000000      ];
       alpha=[   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000      ];
       gamma=[   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000      ];
       beta=[   0.0000000000000000        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.33333333333333331        0.0000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000       0.50000000000000000        0.0000000000000000        0.0000000000000000     
       0.0000000000000000        0.0000000000000000        1.0000000000000000        0.0000000000000000      ];
end
end

function [u,p] = forback1(un,pn,Fun,Fpn,badv,dx,dt,ns)
  global  gfmat 
  gmat = badv/dx*gfmat;
  
  ndofs=size(un,2);
  uu = un';
  pp = pn';
  iter = 0;
  for iter=1:ns
      uu = uu + dt*gmat(1:ndofs,:)*[uu;pp]+dt*Fun';
      pp = pp + dt*gmat(ndofs+1:2*ndofs,:)*[uu;pp]+dt*Fpn';    
  end
  u = uu';
  p = pp'';
end



