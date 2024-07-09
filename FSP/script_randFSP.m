% Script that quantifies the most susceptible faults near Trüllikon.  
% Used to make Figures 3 & S1.
clear;

% Predefine.
dep=1.5;

% Get the fault slip variables.
S=FaultSlipPotential(dep,'deep',1e3,'rand');
str=S.str; dip=S.dip;

% Get average state of stress (and friction).
for i=1:length(S)
    Po(:,:,i)=S(i).Po*(1.12/dep);
    rak(:,:,i)=S(i).rak;
end
PoA=mean(Po,3);
rak=mean(rak,3);
Sh=mean([S.Sh]);
Sv=mean([S.Sv]);
SH=mean([S.SH]);
Pp=mean([S.Pp]);
mu=median([S.mu]);
azi=mean([S.azi]);
if(SH>Sv)
    regime_flag='strike-slip';
else
    regime_flag='normal';
end

% Compute MC parameters.
[STR,DIP]=meshgrid(str,dip);
[Nf,Sf,~]=computeNS(STR(:),DIP(:),Sh,Sv,SH,azi,regime_flag);
%rak=reshape(rak,size(Po));

% Find the most likely slip planes.
[i,j]=find(PoA==min(PoA(:)));
str_o=str(j); dip_o=dip(i); rak_o=rak(sub2ind(size(PoA),i,j))';

% Plot.
figure(1); clf;

% Define an overpressure colormap.
cMap = interp1([0;0.05;0.10;0.25;0.30;0.5;1.0],[0.4 0 0; 1 0 0; 1 0.5 0; 1 0.9 0; 1 1 0; 0 0.95 0; 0 0.1 0],linspace(0,1,256));

% Overpressure fault stereoplot.
%subplot(121);
polarcont(str,dip,PoA,200); hold on;
polarXY(str_o,dip_o,'ow');
polarXY(0:360,30,'-k');
polarXY(0:360,60,'-k');
polarXY(0:360,90,'-k');
polarXY(  0,0:90,'-k');
polarXY( 30,0:90,'-k');
polarXY( 60,0:90,'-k');
polarXY( 90,0:90,'-k');
polarXY(120,0:90,'-k');
polarXY(150,0:90,'-k');
polarXY(180,0:90,'-k');
polarXY(210,0:90,'-k');
polarXY(240,0:90,'-k');
polarXY(270,0:90,'-k');
polarXY(300,0:90,'-k');
polarXY(330,0:90,'-k');
%polarXY(360,0:90,'-k');
colormap(cMap);
clim([max([min(PoA(:)) 0]) max(PoA(:))]);

% Mohr-Coloumb plot.
%subplot(122);
figure(2); clf;
% Circles.
plot([Pp Sh Sv SH],[0 0 0 0],'xk'); hold on;
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
plot(x1,y1,'-k');
plot(x2,y2,'-k');
plot(x3,y3,'-k');
% Friction line.
n=Pp:max(Nf); plot(n,mu*(n-Pp),'-k');
% Overpressure area colors.
[Pc,n,s]=makeGrid([min(y1) max(y1)],[min(x1) max(x1)],Pp,mu,[x2,x3,fliplr(x1)],[y2,y3,fliplr(y1)]); contourf(n,s,Pc,256,'LineColor','none');
% Stress locations of known faults.
%plot(n2,s2,'ok');
% Labels.
xlabel('Normal Stress (MPa)'); ylabel('Shear Stress (MPa)');
colormap(cMap); clim([max([min(PoA(:)) 0]) max(PoA(:))]);
h=colorbar(); ylabel(h, 'Overpressure for Slip (MPa)');
axis equal;
ylim([0 7]);

% Overpressure histogram.
figure(3); clf;
histogram(log10([Po(i(1),j(1),:)]));
xlabel('Overpressure Required to Reactivate the Most Susceptible Fault (log_{10}[MPa])'); ylabel('Counts');








%%%% SUBROUNTINES.

% Get the boundaries of the Mohr-Coloumb circle.
function [x,y]=makeCircle(xl,xr)
  xc=mean([xl xr]);
  radius=abs(xc-xl);
  theta=linspace(0,pi(),200);
  x=radius*cos(theta)+xc;
  y=radius*sin(theta);
  return
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
function [Sn,T,rake]=computeNS(str,dip,Sh,Sv,SH,azi,regime_flag)
  % https://dnicolasespinoza.github.io/node38.html
  
  % Setup problem, depending on stress regime.
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