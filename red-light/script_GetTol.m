clear;

% Predefine some variables.
SMfile='/Users/rschultz/Desktop/CITru/codes/red-light/ShakeMaps/StGallen_M35_20130720/current/products/pgv_withampli.csv'; % StGallen ML 3.5 
%SMfile='/Users/rschultz/Desktop/CITru/codes/red-light/ShakeMaps/Basel_M34_20061208/current/products/pgv_withampli.csv'; % Basel ML 3.4
%SMfile='/Users/rschultz/Desktop/CITru/codes/red-light/ShakeMaps/Basel_M27_20061208/current/products/pgv_withampli.csv'; % Basel ML 2.7
POPfile='/Users/rschultz/Desktop/old/TLP/TLPef/data/pop/LandScan_Global_2018/LSpop.grd';
Nf=1;

% Load in the ShakeMap files.
[latS,lonS,PGV,dPGV]=parseSM(SMfile);
%PGV=PGV/100;

% Get the boundaries of the ShakeMap grid.
latC=[min(latS) max(latS)];
lonC=[min(lonS) max(lonS)];

% Load in population data and grid.
[latG,lonG,POP]=loadPOP(latC,lonC,Nf,POPfile);

% Inpterpolate ShakeMap data onto the population grid.
[LON,LAT]=meshgrid(lonG,latG);
PGV2=griddata(lonS,latS,PGV,LON,LAT,'linear');

% Load in the political boundaries of Switzerland.
S=readgeotable('/Users/rschultz/Desktop/CITru/data/swissboundaries3d_2023-01_2056_5728.shp/swissBOUNDARIES3D_1_4_TLM_LANDESGEBIET.shp');
T = geotable2table(S,["X","Y"]);
[latB,lonB] = projinv(projcrs(2056),T.X{1},T.Y{1});

% Zero off any population not inside Switzerland.
%I = inpolygon(LON(:),LAT(:),lonB,latB);
%POP(~I)=0;

% Compute the chance of being felt.
Pn2=NUISfxn(PGV2,[0 0],2);
Pn3=NUISfxn(PGV2,[0 0],3);
Pn4=NUISfxn(PGV2,[0 0],4);
Pn5=NUISfxn(PGV2,[0 0],5);
Pn6=NUISfxn(PGV2,[0 0],6);

% Compute the number of felt events.
Nn2=sum(POP(:).*Pn2(:));
Nn3=sum(POP(:).*Pn3(:));
Nn4=sum(POP(:).*Pn4(:));
Nn5=sum(POP(:).*Pn5(:));
Nn6=sum(POP(:).*Pn6(:));

% Report those values.
Nn2
Nn3
Nn4

%%

% Plot.
figure(1); clf;
% ShakeMap
subplot(311);
plot3(lonS,latS,PGV); hold on;
plot(lonB,latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title('ShakeMap');
h = colorbar();
% Population data.
subplot(312);
contourf(lonG,latG,log10(POP),'LineColor','none'); hold on;
plot(lonB,latB,'-w','LineWidth',2);
set(gca,'Color','k');
xlabel('Longitude'); ylabel('Latitude'); title(['Population (',sprintf('%0.3g',sum(sum(POP))),')']);
h = colorbar(); colormap(gca,R_colormap('population')); ylabel(h, 'Population (log_{10}[people])'); hold off;
% Interpolated ShakeMap.
subplot(313);
contourf(lonG,latG,PGV2,'LineColor','none'); hold on;
plot(lonB,latB,'-k');
xlabel('Longitude'); ylabel('Latitude'); title('Interpolated ShakeMap');
h = colorbar();


