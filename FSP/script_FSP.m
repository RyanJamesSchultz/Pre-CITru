% Script to make the FSP plots for TRU1-1.
clear;

% Define an overpressure colormap.
cMap = interp1([0;0.05;0.10;0.25;0.30;0.5;1.0],[0.4 0 0; 1 0 0; 1 0.5 0; 1 0.9 0; 1 1 0; 0 0.95 0; 0 0.1 0],linspace(0,1,256));

% Get the fault slip variables.
S=FaultSlipPotential(1.5,'shallow',1,'off');
Sh=S.Sh; Sv=S.Sv; SH=S.SH; azi=S.azi; Pp=S.Pp; regime_flag=S.regime_flag;
mu=S.mu; str=S.str; dip=S.dip;
nf=S.Nf; sf=S.Sf;
Pf=S.Po;
x1=S.x1; y1=S.y1;
x2=S.x2; y2=S.y2;
x3=S.x3; y3=S.y3;

% Get a list of fault orientations of interest.
[str2,dip2,rak2]=getFO();
[n2,s2]=computeNS(str2,dip2,Sh,Sv,SH,azi,regime_flag); 
p2=-(s2-mu*n2+mu*Pp)/mu;



% Plot.
figure(1); clf;

% Overpressure fault stereoplot.
subplot(121);
polarcont(str,dip,Pf,200);
colormap(cMap);
clim([max([min(Pf(:)) 0]) max(Pf(:))]);

% Mohr-Coloumb plot.
subplot(122);
% Circles.
plot([Pp Sh Sv SH],[0 0 0 0],'xk'); hold on;
plot(x1,y1,'-k');
plot(x2,y2,'-k');
plot(x3,y3,'-k');
% Friction line.
n=Pp:max(nf); plot(n,mu*(n-Pp),'-k');
% Overpressure area colors.
[Po,n,s]=makeGrid([min(y1) max(y1)],[min(x1) max(x1)],Pp,mu,[x2,x3,fliplr(x1)],[y2,y3,fliplr(y1)]); contourf(n,s,Po,256,'LineColor','none');
% Stress locations of known faults.
plot(n2,s2,'ok');
% Labels.
xlabel('Normal Stress (MPa)'); ylabel('Shear Stress (MPa)');
colormap(cMap); clim([max([min(Pf(:)) 0]) max(Pf(:))]);
h=colorbar(); ylabel(h, 'Overpressure for Slip (MPa)');
axis equal;




%%%% SUBROUNTINES.

%  Get a list of relevant fault orientations
function [s,d,r]=getFO()
  s=[ 305  217];
  d=[ 059  088];
  r=[-085 -005];
end

% Get the overpressure required for slip inside the Mohr-Coloumb circle.
function [Po,n,s]=makeGrid(Sb,Nb,Pp,f,np,sp)
  n=0.99*min(Nb):0.05:1.01*max(Nb);
  s=0.99*min(Sb):0.05:1.01*max(Sb);
  [N,S]=meshgrid(n,s);
  Po=-(S-f*N+f*Pp)/f;
  I=inpolygon(N,S,np,sp);
  Po(~I)=NaN;
  Po(Po<0)=0;
  return
end

% Get the shear and normal stresses, given a fault dip/strike and the stress field.
function [Sn,T]=computeNS(str,dip,Sh,Sv,SH,azi,regime_flag)
  % https://dnicolasespinoza.github.io/node38.html

  if(strcmpi('reverse',regime_flag))
      S=diag([SH,Sh,Sv]);
      a=azi;
      b=0;
      c=0;
  elseif(strcmpi('strike-slip',regime_flag))
      S=diag([SH,Sv,Sh]);
      a=azi;
      b=0;
      c=90;
  elseif(strcmpi('normal',regime_flag))
      S=diag([Sv,SH,Sh]);
      a=azi-90;
      b=90;
      c=0;
  end
  R=[cosd(a)*cosd(b),                         sind(a)*cosd(b),                         -sind(b);
     cosd(a)*sind(b)*sind(c)-sind(a)*cosd(c), sind(a)*sind(b)*sind(c)+cosd(a)*cosd(c), cosd(b)*sind(c);
     cosd(a)*sind(b)*cosd(c)+sind(a)*sind(c), sind(a)*sind(b)*cosd(c)-cosd(a)*sind(c), cosd(b)*cosd(c)];
  Sg=R'*S*R;
  
  Sn=zeros(size(str));
  T=Sn;
  
  for i=1:length(str)

      nn=[-sind(str(i))*sind(dip(i)); cosd(str(i))*sind(dip(i)); -cosd(dip(i))];
      ns=[ cosd(str(i));              sind(str(i));               0           ];
      nd=[-sind(str(i))*cosd(dip(i)); cosd(str(i))*cosd(dip(i));  sind(dip(i))];
  
      t=Sg*nn;
      Sn(i)=dot(t,nn);

      Td=dot(t,nd);
      Ts=dot(t,ns);
    
      T(i)=sqrt(Ts^2+Td^2);
      rake=atan2d(Td,Ts);
  end
  
  return
end



