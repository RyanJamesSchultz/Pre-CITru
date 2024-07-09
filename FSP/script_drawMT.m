% Script to plot earthquake focal mechanisms.
% Assists with making Figures 3 & S1.
clear;

% Get the strike-dip angles of the faults.
str=[168.00 348.00  14 194];
dip=[ 52.00  52.00  90  90];
rak=[-90.45 -90.45   0   0];
mag=[  3      3      3   3];

% Plot moment tensors.
figure(4); clf;
x=4*(1:length(mag));
for i=1:length(mag)
    mt=getMT(str(i),dip(i),rak(i));
    focalmech(mt,x(i),0,mag(i));
end
axis equal;
set(gca,'Visible','off');



%%%% SUBROUNTINE.

% Get the moment tensor from strike/dip/rake.
function [M]=getMT(strike,dip,rake)
  
  % Convert to Cartesian tensor elements (Aki & Richards (2002), p 112-113).
  Mxx=-( sind(dip)*cosd(rake)*sind(2*strike)  +     sind(2*dip)*sind(rake)*(sind(strike)^2) );
  Mxy=+( sind(dip)*cosd(rake)*cosd(2*strike)  + 0.5*sind(2*dip)*sind(rake)*sind(2*strike)   );
  %Myx=Mxy;
  Mxz=-( cosd(dip)*cosd(rake)*cosd(strike)    +     cosd(2*dip)*sind(rake)*sind(strike)     );
  %Mzx=Mxz;
  Myy=+( sind(dip)*cosd(rake)*sind(2*strike)  -     sind(2*dip)*sind(rake)*(cosd(strike)^2) );
  Myz=-( cosd(dip)*cosd(rake)*sind(strike)    -     cosd(2*dip)*sind(rake)*cosd(strike)     );
  %Mzy=Myz;
  Mzz=sind(2*dip)*sind(rake);

  % GMT CMT formatting
  Mrr=Mzz;
  Mrdelta=Mxz;
  Mrphi=-Myz;
  Mdeltadelta=Mxx;
  Mdeltaphi=-Mxy;
  Mphiphi=Myy;

  % Stuff together.
  % moment tensor (mrr, mtt, mff, mrt, mrf, mtf).
  M=[Mrr,Mdeltadelta,Mphiphi,Mrdelta,Mrphi,Mdeltaphi];
  return
end

