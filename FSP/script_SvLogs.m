clear;

% Pulled Sv data from Figure 5-20 of NAGRA's TRU1-1 report, Dossier VI, Page 64 (p.72).
z=[100 200 300 400 500 600 700 800 900 1000 1100 1200 1300];
Sv=[16.8318 33.9095 50.7361 67.7488 85.0563 103.9125 122.9567 141.5391 159.8102 178.298 197.3812 216.9817 235.9056]*35/255.385;

% Plot results.
figure(1); clf;
subplot(121);
plot(Sv,z);
xlabel('Vertical Stress, S_V (MPa)'); ylabel('Depth (m)');
set(gca, 'YDir','reverse');
subplot(122);
plot(Sv./(z*1000),z);
xlabel('Vertical Stress Gradient, S_V (MPa/km)'); ylabel('Depth (m)');
set(gca, 'YDir','reverse');

% Report values.
%mean()

