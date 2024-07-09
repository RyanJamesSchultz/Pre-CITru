function S=FaultSlipPotential(z,depth_flag,N,rand_flag)
  % Function to get the amount of overpressure required to initiate fault
  % slip, for a given stress field.
  % 
  % Written by Ryan Schultz.
  % 
  
  % Predfeine some values.
  S=struct('Sh',[],'Sv',[],'SH',[],'azi',[],'Pp',[],'regime_flag',[],'mu',[],'str',[],'dip',[],'rak',[],'Nf',[],'Sf',[],'Po',[],'x1',[],'y1',[],'x2',[],'y2',[],'x3',[],'y3',[]);
  str=0:360; dip=0:90;
  [STR,DIP]=meshgrid(str,dip);
  
  % Loop over all of the iterations
  for i=1:N
      
      % Output for percent done.
      100*i/N
      
      % Define the stress tensor and friction value.
      [Sh,Sv,SH,azi,Pp,regime_flag]=getS(z,depth_flag,rand_flag);
      mu=getMU(Sh,Sv,SH,Pp,regime_flag);
      
      % Get the overpressure values for all possible fault orientations.
      [Nf,Sf,rak]=computeNS(STR(:),DIP(:),Sh,Sv,SH,azi,regime_flag); 
      Po=-(Sf-mu*Nf+mu*Pp)/mu;
      Po=reshape(Po,size(STR));
      rak=reshape(rak,size(STR));
      
      % Get Mohr-Coloumb cirlces.
      if(strcmpi('reverse',regime_flag))
          [x1,y1]=makeCircle(Sv,SH); 
          [x2,y2]=makeCircle(Sv,Sh); 
          [x3,y3]=makeCircle(Sh,SH); 
      elseif(strcmpi('strike-slip',regime_flag))
          [x1,y1]=makeCircle(Sh,SH); 
          [x2,y2]=makeCircle(Sv,SH); 
          [x3,y3]=makeCircle(Sh,Sv); 
      elseif(strcmpi('normal',regime_flag))
          [x1,y1]=makeCircle(Sh,Sv); 
          [x2,y2]=makeCircle(Sv,SH); 
          [x3,y3]=makeCircle(Sh,SH); 
      end
      
      % Save everything in a structure.
      S(i).Sh=Sh; S(i).Sv=Sv; S(i).SH=SH; S(i).azi=azi; S(i).Pp=Pp; S(i).regime_flag=regime_flag;
      S(i).mu=mu; S(i).str=str; S(i).dip=dip; S(i).rak=rak; S(i).Sf=Sf; S(i).Nf=Nf; S(i).Po=Po;
      S(i).x1=x1; S(i).y1=y1; S(i).x2=x2; S(i).y2=y2; S(i).x3=x3; S(i).y3=y3; 
      
  end

  return
end




%%%% SUBROUNTINES.

% Get stress state information from prior study models.
function [Sh,Sv,SH,azi,Pp,regime_flag]=getS(z,dep_flag,rand_flag)
  % Diehl, Madritsch, Schnellmann, Spillmann, Brockmann, & Wiemer (2023). Seismotectonic evidence for present-day transtensional reactivation of the slowly deforming Hegau-Bodensee Graben in the northern foreland of the Central Alps. Tectonophysics, 846, 229659, doi: 10.1016/j.tecto.2022.229659.
  % Gonus, Bailey, Desroches, & Garrard (2021). TRU - Dossier VI Wireline Logging and Microhydraulic Fracturing. NAGRA Report.
  % Garrard, Gonus, Desroches, & Bailey (2021). BUL - Dossier VI Wireline Logging and Microhydraulic Fracturing. NAGRA Report.

  % Get the state of stress.
  Sh=19.40*z-00.0; % From TRU1-1 & BUL1-1 reports.
  Sv=24.84*z-00.0; % From integration of density to basement (TRU1-1 report).
  Pp= 9.81*z-00.0; % Hydrostatic pore-pressure.
  
  % Perturb values, if flagged to.
  if(strcmpi(rand_flag,'rand'))
      Sh=Sh+normrnd(0,3.1,[1 1]);
      Sv=Sv+normrnd(0,0.1,[1 1]);
      Pp=Pp+normrnd(0,1.0,[1 1]);
  end
  
  % Stress regime changes with depth.
  if(strcmpi(dep_flag,'shallow'))
      azi=78+90;
      R  =mean([0.27 0.19]); % Shape ratio from EQ moment tensors.
      phi=mean([0.73 0.81]);
      Aphi=phi;
      if(strcmpi(rand_flag,'rand'))
          azi =azi +normrnd(0,4.0,[1 1]);
          Aphi=Aphi+normrnd(0,0.2,[1 1]);
      end
  elseif(strcmpi(dep_flag,'deep'))
      azi=247-90;
      R  =0.13; % Shape ratio from EQ moment tensors.
      phi=0.87;
      Aphi=2-phi;
      if(strcmpi(rand_flag,'rand'))
          azi =azi +normrnd(0,3.0,[1 1]);
          Aphi=Aphi+normrnd(0,0.2,[1 1]);
      end
  end
  
  % Get the value of SHmax and the stress regime.
  if(Aphi<1)
      phi=Aphi;
      R=1-phi;
      SH=R*(Sh-Sv)+Sv;
      regime_flag='normal';
  elseif(Aphi>2)
      phi=Aphi-2;
      R=1-phi;
      SH=(R*Sv-Sh)/(R-1);
      regime_flag='reverse';
  else
      phi=2-Aphi;
      R=1-phi;
      SH=(R*Sh-Sv)/(R-1);
      regime_flag='strike-slip';
  end
  
  return
end


% Get the equilibrium friction value, given the stress state.
function [mu]=getMU(Sh,Sv,SH,Pp,regime_flag)
  
  % Sus out stress tensor ordering, from stress regime.
  if(strcmpi(regime_flag,'reverse'))
      S3=Sv; S2=Sh; S1=SH;
  elseif(strcmpi(regime_flag,'strike-slip'))
      S3=Sh; S2=Sv; S1=SH;
  elseif(strcmpi(regime_flag,'normal'))
      S3=Sh; S2=SH; S1=Sv;
  end
  
  % Compute equlibrium friction value.
  nc=mean([S1,S3]);
  r=S1-nc;
  a=(8*Pp*nc-4*nc^2+4*r^2-4*Pp^2);
  c=(4*r^2);
  mu=sqrt(-c)/sqrt(a);

  return
end


% Get the shear and normal stresses, given a fault dip/strike and the stress field.
function [Sn,T,rake]=computeNS(str,dip,Sh,Sv,SH,azi,regime_flag)
  % https://dnicolasespinoza.github.io/node38.html
  
  % Setup problem, depending on stress regime.
  if(strcmpi(regime_flag,'reverse'))
      S=diag([SH,Sh,Sv]);
      a=azi;
      b=0;
      c=0;
  elseif(strcmpi(regime_flag,'strike-slip'))
      S=diag([SH,Sv,Sh]);
      a=azi;
      b=0;
      c=90;
  elseif(strcmpi(regime_flag,'normal'))
      S=diag([Sv,SH,Sh]);
      a=azi-90;
      b=90;
      c=0;
  end
  R=[cosd(a)*cosd(b),                         sind(a)*cosd(b),                         -sind(b);
     cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c), sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c), cosd(b)*sind(c);
     cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c), sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c), cosd(b)*cosd(c)];
  Sg=R'*S*R;
  
  % Transform into the the fault's coordinates.
  Sn=zeros(size(str));
  T=Sn; rake=Sn;
  for i=1:length(str)

      nn=[-sind(str(i))*sind(dip(i)); cosd(str(i))*sind(dip(i)); -cosd(dip(i))];
      ns=[ cosd(str(i));              sind(str(i));               0           ];
      nd=[-sind(str(i))*cosd(dip(i)); cosd(str(i))*cosd(dip(i));  sind(dip(i))];
  
      t=Sg*nn;
      Sn(i)=dot(t,nn);

      Td=dot(t,nd);
      Ts=dot(t,ns);
    
      T(i)=sqrt(Ts^2+Td^2);
      rake(i)=atan2d(Td,Ts);
  end
  
  return
end


% Get the boundaries of the Mohr-Coloumb circle.
function [x,y]=makeCircle(xl,xr)
  xc=mean([xl xr]);
  radius=abs(xc-xl);
  theta=linspace(0,pi(),200);
  x=radius*cos(theta)+xc;
  y=radius*sin(theta);
  return
end

