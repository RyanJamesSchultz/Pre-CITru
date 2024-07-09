% Script that quantifies the slip potential of faults near Trüllikon.  
% Used to make Figure 4.
clear;

% Predefine.
dep=1.5;

% Load in all of the data.
[vID1,vP1,tIDs1]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_Neuhausen-1.ts');
[vID2,vP2,tIDs2]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_Neuhausen-2.ts');
[vID3,vP3,tIDs3]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_Neuhausen-3.ts');
[vID4,vP4,tIDs4]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_Stadel_Irchel_1.ts');
[vID5,vP5,tIDs5]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_Stadel_Irchel_7.ts');
[vID6,vP6,tIDs6]=parseTS('/Users/rschultz/Desktop/Pre-CITru/data/faults_horizons_CITRU/faults/selected/file_ST_Fault_ZNO1.ts');
TRU1 = [695372.648,277548.076,475.07; 695372.648,277548.076,475.07-(1084+35)];

% Get the vertex information.
[tM1,tN1,str1,dip1]=getTriangleData(tIDs1,vID1,vP1); % NF1
[tM2,tN2,str2,dip2]=getTriangleData(tIDs2,vID2,vP2); % NF2
[tM3,tN3,str3,dip3]=getTriangleData(tIDs3,vID3,vP3); % NF3
[tM4,tN4,str4,dip4]=getTriangleData(tIDs4,vID4,vP4); % BIH-FZ1
[tM5,tN5,str5,dip5]=getTriangleData(tIDs5,vID5,vP5); % BIH-FZ2
[tM6,tN6,str6,dip6]=getTriangleData(tIDs6,vID6,vP6); % OF

% Get the fault slip variables.
S=FaultSlipPotential(dep,'shallow',1e3,'rand');
str=S.str; dip=S.dip;

% Get average state of stress (and friction).
for i=1:length(S)
    Po(:,:,i)=S(i).Po;
    rak(:,:,i)=S(i).rak;
end
Po=mean(Po,3);
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

% Compute the overpressure for each fault vertex.
op1=interp2(str,dip,Po,str1,dip1,'linear'); % NF1
op2=interp2(str,dip,Po,str2,dip2,'linear'); % NF2
op3=interp2(str,dip,Po,str3,dip3,'linear'); % NF3
op4=interp2(str,dip,Po,str4,dip4,'linear'); % BIH-FZ1
op5=interp2(str,dip,Po,str5,dip5,'linear'); % BIH-FZ2
op6=interp2(str,dip,Po,str6,dip6,'linear'); % OF

% Change reference frame to use the TRU1-1 bottom-hole as an origin.
tM1=tM1-TRU1(2,:);  vP1=vP1-TRU1(2,:);
tM2=tM2-TRU1(2,:);  vP2=vP2-TRU1(2,:);
tM3=tM3-TRU1(2,:);  vP3=vP3-TRU1(2,:);
tM4=tM4-TRU1(2,:);  vP4=vP4-TRU1(2,:);
tM5=tM5-TRU1(2,:);  vP5=vP5-TRU1(2,:);
tM6=tM6-TRU1(2,:);  vP6=vP6-TRU1(2,:);
TRU1=TRU1-TRU1(2,:);

% Change to units of km, from m.
tM1=tM1/1e3;  vP1=vP1/1e3;
tM2=tM2/1e3;  vP2=vP2/1e3;
tM3=tM3/1e3;  vP3=vP3/1e3;
tM4=tM4/1e3;  vP4=vP4/1e3;
tM5=tM5/1e3;  vP5=vP5/1e3;
tM6=tM6/1e3;  vP6=vP6/1e3;
TRU1=TRU1/1e3;

% Make a second coordinate system oriented along the average strike
% direction of the Neuhausen Fault.
theta=mean([str1;str2;str3])-15;
%theta=-45;
tM1p=Rotate(tM1,theta);  vP1p=Rotate(vP1,theta);
tM2p=Rotate(tM2,theta);  vP2p=Rotate(vP2,theta);
tM3p=Rotate(tM3,theta);  vP3p=Rotate(vP3,theta);
tM4p=Rotate(tM4,theta);  vP4p=Rotate(vP4,theta);
tM5p=Rotate(tM5,theta);  vP5p=Rotate(vP5,theta);
tM6p=Rotate(tM6,theta);  vP6p=Rotate(vP6,theta);
TRU1p=Rotate(TRU1,theta);

% Compute the distance from TRU1-1 bottomhole to each vertex.
R1=vecnorm(tM1')';
R2=vecnorm(tM2')';
R3=vecnorm(tM3')';
R4=vecnorm(tM4')';
R5=vecnorm(tM5')';
R6=vecnorm(tM6')';

% Closest point on the Neuhausen Fault to TRU1-1.
[~,Ic]=min(R3);
for i=1:length(S)
    opC(i)=interp2(S(i).str,S(i).dip,S(i).Po,str3(Ic),dip3(Ic),'linear')*(1.12/dep);
end

% Most susceptible point on the Neuhausen Fault.
temp=op3; %temp(R3>2)=Inf; %temp(R3>2)=Inf;
[~,Is]=min(temp);
for i=1:length(S)
    opS(i)=interp2(S(i).str,S(i).dip,S(i).Po,str3(Is),dip3(Is),'linear')*(1.12/dep);
end

% Close/susceptible point on the N-S fault.
[~,Ie]=min(op6);
for i=1:length(S)
    opE(i)=interp2(S(i).str,S(i).dip,S(i).Po,str6(Ie),dip6(Ie),'linear')*(1.12/dep);
end

% Temp change for finding most susceptible part to slip.
%I=tM3p(:,3)<0;
%R3=R3(I);
%tM3p=tM3p(I,:);
%tM3=tM3(I,:);
%op3=op3(I);

% Define an overpressure colormap.
cMap = interp1([0;0.05;0.10;0.25;0.30;0.5;1.0],[0.4 0 0; 1 0 0; 1 0.5 0; 1 0.9 0; 1 1 0; 0 0.95 0; 0 0.1 0],linspace(0,1,256));

% Plot Figure 4.
figure(4); clf;
% Map.
subplot(231); hold on;
scatter(tM1(:,1),tM1(:,2),20,op1,'o','filled');
scatter(tM2(:,1),tM2(:,2),20,op2,'o','filled');
scatter(tM3(:,1),tM3(:,2),20,op3,'o','filled');
scatter(tM4(:,1),tM4(:,2),20,op4,'o','filled');
scatter(tM5(:,1),tM5(:,2),20,op5,'o','filled');
scatter(tM6(:,1),tM6(:,2),20,op6,'o','filled');
plot(TRU1(:,1),TRU1(:,2),'-ok','MarkerFaceColor','k');
xlabel('Easting (km)'); ylabel('Northing (km)'); zlabel('Depth (km)');

% Stereoplot.
subplot(234); hold on;
polarcont(str,dip,Po,200);
polarXY(str1,dip1,'.','MarkerEdgeColor','b');
polarXY(str2,dip2,'.','MarkerEdgeColor','b');
polarXY(str3,dip3,'.','MarkerEdgeColor','b');
polarXY(str4,dip4,'.','MarkerEdgeColor','#4B0082');
polarXY(str5,dip5,'.','MarkerEdgeColor','#4B0082');
polarXY(str6,dip6,'.','MarkerEdgeColor','k');
colormap(cMap);
clim([max([min(Po(:)) 0]) max(Po(:))]);

% Fault-parallel depth profile.
subplot(2,3,[2 3 5 6]); hold on;
scatter(tM1p(:,1),tM1p(:,3),20,op1,'o','filled');
scatter(tM2p(:,1),tM2p(:,3),20,op2,'o','filled');
scatter(tM3p(:,1),tM3p(:,3),20,op3,'o','filled');
scatter(tM4p(:,1),tM4p(:,3),20,op4,'o','filled');
scatter(tM5p(:,1),tM5p(:,3),20,op5,'o','filled');
scatter(tM6p(:,1),tM6p(:,3),20,op6,'o','filled');
plot(TRU1p(:,1),TRU1p(:,3),'-ok','MarkerFaceColor','k');
xlabel('Fault Parallel Direction (km)'); ylabel('Elevation (km)');
colormap(cMap);
%%
% 3D Plot (in N-E-Z coords).
figure(1); clf;
plot3(TRU1(:,1),TRU1(:,2),TRU1(:,3),'-ok','MarkerFaceColor','k'); hold on;
% scatter3(vP1(:,1),vP1(:,2),vP1(:,3),'o','MarkerEdgeColor','b');
% scatter3(tM1(:,1),tM1(:,2),tM1(:,3),20,op1,'o','filled');
% scatter3(vP2(:,1),vP2(:,2),vP2(:,3),'o','MarkerEdgeColor','b');
% scatter3(tM2(:,1),tM2(:,2),tM2(:,3),20,op2,'o','filled');
% scatter3(vP3(:,1),vP3(:,2),vP3(:,3),'o','MarkerEdgeColor','b');
% scatter3(tM3(:,1),tM3(:,2),tM3(:,3),20,op3,'o','filled');
% scatter3(tM3(Ic,1),tM3(Ic,2),tM3(Ic,3),'ob','filled');
% scatter3(tM3(Is,1),tM3(Is,2),tM3(Is,3),'ob','filled');
% scatter3(vP4(:,1),vP4(:,2),vP4(:,3),'o','MarkerEdgeColor','#4B0082');
% scatter3(tM4(:,1),tM4(:,2),tM4(:,3),20,op4,'o','filled');
% scatter3(vP5(:,1),vP5(:,2),vP5(:,3),'o','MarkerEdgeColor','#4B0082');
% scatter3(tM5(:,1),tM5(:,2),tM5(:,3),20,op5,'o','filled');
scatter3(vP6(:,1),vP6(:,2),vP6(:,3),'o','MarkerEdgeColor','k');
scatter3(tM6(:,1),tM6(:,2),tM6(:,3),20,op6,'o','filled');
scatter3(tM6(Ie,1),tM6(Ie,2),tM6(Ie,3),'ob','filled');
xlabel('Easting (km)'); ylabel('Northing (km)'); zlabel('Depth (km)');
colormap(cMap);

% Stress plots.
figure(2); clf;

% Overpressure fault stereoplot.
subplot(121); hold on;
polarcont(str,dip,Po,200);
polarXY(str1,dip1,'o','MarkerEdgeColor','b');
polarXY(str2,dip2,'o','MarkerEdgeColor','b');
polarXY(str3,dip3,'o','MarkerEdgeColor','b');
polarXY(str4,dip4,'o','MarkerEdgeColor','#4B0082');
polarXY(str5,dip5,'o','MarkerEdgeColor','#4B0082');
polarXY(str6,dip6,'o','MarkerEdgeColor','k');
colormap(cMap);
clim([max([min(Po(:)) 0]) max(Po(:))]);

% Mohr-Coloumb plot.
subplot(122); hold on;
% Circles.
plot([Pp Sh Sv SH],[0 0 0 0],'xk');
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
colormap(cMap); clim([max([min(Po(:)) 0]) max(Po(:))]);
h=colorbar(); ylabel(h, 'Overpressure for Slip (MPa)');
axis equal;

% R vs OP plot.
figure(3); clf; hold on;
plot(R1,op1,'o','MarkerEdgeColor','b');
plot(R2,op2,'o','MarkerEdgeColor','b');
plot(R3,op3,'o','MarkerEdgeColor','b');
plot(R3(Ic),op3(Ic),'x','MarkerEdgeColor','b');
plot(R3(Is),op3(Is),'x','MarkerEdgeColor','b');
plot(R4,op4,'o','MarkerEdgeColor','#4B0082');
plot(R5,op5,'o','MarkerEdgeColor','#4B0082');
plot(R6,op6,'o','MarkerEdgeColor','k');
plot(R6(Ie),op6(Ie),'x','MarkerEdgeColor','k');
xlabel('Distance (km)'); ylabel('Overpressure (MPa)')

% Histogram of overpressures at the closest and most susceptbile points of
% the Neuhausen Fault.
figure(5); clf;
subplot(211);
histogram(opC);
set(gca,'xscale','log');
xlabel('Overpressure (MPa)'); ylabel('Counts');
subplot(212);
histogram(opS);
set(gca,'xscale','log');
xlabel('Overpressure (MPa)'); ylabel('Counts');

% Reporting average values of overpressure gradient on the faults.
disp('Average FSP (NF & BIH-FZ)')
mean([op1;op2;op3])/dep
nanmean([op4;op5])/dep

% Report percentage of cases encountering IS at the Neuhausen Fault.
disp('FSP Chances')
100*sum(opC<=0.0475)/length(opC)
100*sum(opS<=0.0260)/length(opS)
100*sum(opE<=0.0445)/length(opE)

%%%% SUBROUNTINES.

% Get the midpoint of a vertex triangle.
function [tM,tN,str,dip]=getTriangleData(tID,vID,vP)
  
  n=size(tID,1);
  tM=zeros([n 3]);
  tN=zeros([n 3]);
  str=zeros([n 1]);
  dip=zeros([n 1]);
  
  for i=1:n
      
      % Get the three vertexes for this triangle.
      tid=tID(i,:);
      I=find(any((vID==tid)'));
      vp=vP(I,:); % N-S,E-W,U-D vector.
      
      % Find the mid-point of the three vertexes.
      tM(i,:)=mean(vp);
      
      % Find the normal vector of the three vertexes.
      tN(i,:)=cross(vp(2,:)-vp(1,:),vp(3,:)-vp(1,:));
      tN(i,:)=tN(i,:)/norm(tN(i,:));
      if(tN(i,3)<0)
          tN(i,:)=-tN(i,:);
      end
      
      % Find the strike and dip.
      dip(i)=acosd(tN(i,3));
      str(i)=atan2d(tN(i,1),tN(i,2))-90;
      
  end
  str(str<0)=str(str<0)+360;
  
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


% Rotate to a new coordinate system.
function [Xp]=Rotate(X,theta)
  Xp=zeros(size(X));
  R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)];
  for i=1:size(X,1)
      Xp(i,1:2)=X(i,1:2)*R;
      Xp(i,3)=X(i,3);
  end
  return
end

