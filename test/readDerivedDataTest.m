clear; close all; clc;

%% Initializations
constellations = {'G' 'S' 'R' 'J' 'C' 'E' '' '' '-'};
ID_GPS  = find(strcmp(constellations, 'G'));
ID_SBAS = find(strcmp(constellations, 'S'));
ID_GLO  = find(strcmp(constellations, 'R'));
ID_QZSS = find(strcmp(constellations, 'J'));
ID_BDS  = find(strcmp(constellations, 'C'));
ID_GAL  = find(strcmp(constellations, 'E'));

fL1 = 1.57542e9;
fL5 = 1.17645e9;

%% Configuration
% Measurement campaign name
campaignName = '2020-05-21-US-MTV-1';
% Phone name
phoneName = 'Pixel4XL';
% Constellation to visualize
constel = ID_GPS;
% SV ID (prn) to visualize
svId = 30;
% Frequency band to visualize
fBand = fL5;

%% Dataset
% Dataset path
datasetsPath = './data/training/datasets/';
testFileName = [phoneName '_GnssLog.derived'];
dirName = [datasetsPath campaignName '/'];

%% Read data
Derived = readDerivedData(dirName, testFileName);

%% Testing plots
figure;
histogram(Derived.MEAS.ConstellationType);
xticklabels(constellations);
xlabel('Constellation')

constelData = Derived.MEAS(Derived.MEAS.ConstellationType == constel, :);

figure;
histogram(constelData.Svid);
xlabel('SV ID (GPS)')

satData = constelData((constelData.Svid == svId & constelData.CarrierFrequencyHz == fBand), :);
figure;
hold on;
plot(satData.TimeNanos/1e9, satData.RawPrM);
plot(satData.TimeNanos/1e9, satData.SmPrM);

xlabel('Time (s)'); ylabel('Pseudorange (m)')
legend('Raw', 'Smoothed');

figure;
plot(satData.TimeNanos/1e9, satData.ElDeg);
xlabel('Time (s)'); ylabel('Elevation (deg)')

figure;
hold on;
plot(satData.TimeNanos/1e9, satData.RawPrErrorM);
plot(satData.TimeNanos/1e9, satData.SmPrErrorM);
xlabel('Time (s)'); ylabel('Pseudorange error (m)')
legend('Raw', 'Smoothed');

pctl = [50 95];
pctlRaw = prctile(abs(satData.RawPrErrorM), pctl);
pctlSm = prctile(abs(satData.SmPrErrorM), pctl);
figure;
hold on;
ecdf(abs(satData.RawPrErrorM));
ecdf(abs(satData.SmPrErrorM));
for k = 1:numel(pctl)
    plot([1;1]*pctlRaw(k), [0;1]*pctl(k)/100, '--b')
    plot([0;1]*pctlRaw(k), [1;1]*pctl(k)/100, '--b')
    plot([1;1]*pctlSm(k), [0;1]*pctl(k)/100, '--r')
    plot([0;1]*pctlSm(k), [1;1]*pctl(k)/100, '--r')
end
xlabel('Pseudorange error (m)')
legend('Raw', 'Smoothed');




