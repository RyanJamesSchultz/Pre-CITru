clear;

dS=3e6;  % Pa
L=20000; w=10000;  % NF (m).
L= 2100; w=  600;  % NS (m).

Mo=dS*w^2*L*pi()/2; % N m.
Mw=(log10(Mo)-9.1)/1.5; 
%Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).