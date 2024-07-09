% Script that quantifies the red-light for Trüllikon.
% Used to make Figures 7, S3, & S4.
clear;

% Predefine some variables.
Rfile='RiskCurve.csv';
Nfile='NuisCurve.csv';
%b=1.0;
b=0.7; % for FS3.
%b=1.3; % for FS3.
Mmin=0.0;
Mmax=6.6; % NF.
%Mmax=4.5; % NS. for FS2.

% Predefine the risk tolerances.
TdG2=2.1622;
TnG1=1e4;
TnG2=1.4102e+05;
TnB1=2.5622e+04;
TnB2=1.7566e+05;
TnUK=5478*2.4;
TnAB=3.5e4*4.0;

% Report the mid points.
geomean([TnB1 TnB2])
geomean([TnG1 TnG2])
geomean([geomean([TnB1 TnB2]) geomean([TnG1 TnG2])])

% Get the risk data.
[  Rd,Rc,Ri,Rf,Re,Mr]=parseRISK(Rfile);
[~,N2,N3,N4,N5,N6,Mn]=parseNUIS(Nfile);

% Convert from number of people who felt into number of homes exposed.
%N2=N2/4; N3=N3/4; N4=N4/4; N5=N5/4; N6=N6/4;

% Interpolate the scenario curves.
M2=Mmin:0.1:max(Mr);
Re2=10.^interp1(Mr,log10(Re),M2,'pchip','extrap'); Re2(Re2<1)=0; Re2(M2<3.5)=0;
Rf2=10.^interp1(Mr,log10(Rf),M2,'pchip','extrap'); Rf2(Rf2<1)=0;
Rn2=10.^interp1(Mn,log10(N3),M2,'pchip','extrap'); Rn2(Rn2<1)=0;

% Interpolate the scenario curves.
Mo=Mmin:0.01:Mmax;
Re3=10.^interp1(Mr,log10(Re),Mo,'pchip','extrap'); Re3(Re3<1)=0; Re3(Mo<3.5)=0;
Rf3=10.^interp1(Mr,log10(Rf),Mo,'pchip','extrap'); Rf3(Rf3<1)=0;
Rn3=10.^interp1(Mn,log10(N3),Mo,'pchip','extrap'); Rn3(Rn3<1)=0;

% Get the GR-MFD weights.
Wobs=GR_MFD(Mo,[Mmin Mmax],[1 b],'norm');

% Compute the NRBE risks.
for i=1:length(Mo)
    
    Wi=Wobs(i:end);
    Wi=Wi/sum(Wi);
    
    Re_n(i)=getRISKsamp(Mo(i:end),Wi,Re3(i:end));
    Rf_n(i)=getRISKsamp(Mo(i:end),Wi,Rf3(i:end));
    Rn_n(i)=getRISKsamp(Mo(i:end),Wi,Rn3(i:end));
end
Re_n(Re_n<1)=0;
Rf_n(Rf_n<1)=0;
Rn_n(Rn_n<1)=0;

% Plot.
figure(1); clf;

% Plot the scenario risk curves.
subplot(131);
semilogy(M2,Rn2,'-b','DisplayName','Scenario Earthquake'); hold on;
semilogy(Mo,Rn_n,'--b','DisplayName','Next Largest Earthquake');
semilogy(xlim(),TnB1*[1 1],'--k','DisplayName','Basel M_L 2.7');
semilogy(xlim(),TnB2*[1 1],'--k', 'DisplayName','Basel M_L 3.4');
%semilogy(xlim(),TnG1*[1 1],'--b','DisplayName','StGallen ML 2.1');
semilogy(xlim(),TnG2*[1 1],'-k', 'DisplayName','StGallen M_L 3.5');
%semilogy(xlim(),TnUK*[1 1],':k', 'LineWidth',2,'DisplayName','UK tolerance estimate');
%semilogy(xlim(),TnAB*[1 1],':k', 'LineWidth',2,'DisplayName','AB tolerance estimate');
ylabel('Number of People who of Felt (at CDI 3)');
xlabel('Earthquake Magnitude (Mw)');
legend('Location','northwest');
title('Nuisance Risks');
xlim([Mmin max(Mr)]);
grid on;
subplot(132);
semilogy(M2,Re2/1e6,'-r','DisplayName','Scenario Earthquake'); hold on;
semilogy(Mo,Re_n/1e6,'--r','DisplayName','Next Largest Earthquake');
semilogy(xlim(),TdG2*[1 1],'-k','DisplayName','StGallen M_L 3.5');
ylabel('Economic Loss (M CHF)');
xlabel('Earthquake Magnitude (Mw)');
legend('Location','northwest');
title('Damage Losses');
xlim([Mmin max(Mr)]);
grid on;
subplot(133);
semilogy(M2,Rf2,'-g','DisplayName','Scenario Earthquake'); hold on;
semilogy(Mo,Rf_n,'--g','DisplayName','Next Largest Earthquake');
ylabel('Number of Fatalities');
xlabel('Earthquake Magnitude (Mw)');
legend('Location','northwest');
title('Fatality Risks');
xlim([Mmin max(Mr)]);
grid on;






%%% SUBROUTINES.
function R2=getRISKsamp(m,w,R)
    
    % Trivial computation.
    if(length(w)==1)
        R2=R;
        return
    end

    % Mean.
    R2=sum(w.*R); % Mean.
    
    % Median.
    %w=cumsum(w);
    %R2=interp1(w,R,0.5,'linear'); 

end