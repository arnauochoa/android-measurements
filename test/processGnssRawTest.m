clearvars %-except gnssRaw gnssAnalysis accRaw gyrRaw magRaw sensorAnalysis obsRef obsTypeRef obsExtraRef; 
close all; clc;

%% Constants
LIGHTSPEED = 299792458;

%% Configuration
% Measurement campaign name
campaignName = '2021-04-29-US-MTV-1'; % '2020-06-05-US-MTV-1'; '2020-08-06-US-MTV-2';
% Phone name
phoneName = 'Pixel4';
% Filter flag: set to 1 to apply filters to measurements
filter = 1;

%% Dataset
% Dataset path
datasetsPath = [workspacePath 'data/sdc-data/train/'];
rawFileName = [phoneName '_GnssLog.txt'];
% derivedFileName = [phoneName '_GnssLog.derived'];
dirName = [datasetsPath campaignName '/' phoneName filesep];
% rinexFilePath = [dirName phoneName '_GnssLog.20o'];

global figsPath 
figsPath = ['./figs/' campaignName '/' phoneName '/'];
if ~exist(figsPath, 'dir')
    mkdir(figsPath);
end

% %% Read data with Thomas' code
% if ~exist('obsRef', 'var')
%     [obsRef, obsTypeRef, obsExtraRef] = read_AndroidRawMeasRoBot(rawFileName, dirName, 0);
% end

%% Read derived data (from GNSS Analysis)
% Derived = readDerivedData(dirName, derivedFileName);

%% Read RINEX
% disp('Reading RINEX file...');
% [obs_rinex,obs_rinex_type] = rinex_v3_obs_parser(rinexFilePath);

%% Read Android Raw Log
if ~exist('gnssRaw', 'var')
    disp('Reading Android Raw Log...')
    [gnssRaw, gnssAnalysis, accRaw, gyrRaw, magRaw, sensorAnalysis] = readGnssLog(dirName,rawFileName);
    disp('Android Raw Log read.')
end

%% Process Android GNSS raw measurements
[obs, obsType] = processGnssRaw(gnssRaw, filter);

%% Plot IMU + Mag
% figure; plot(accRaw.utcTimeMillis, [accRaw.UncalAccelXMps2 accRaw.UncalAccelYMps2 accRaw.UncalAccelZMps2], '.'); 
% xlabel('Time (ms)'); ylabel('Acceleration (m/s²)');
% legend('X', 'Y', 'Z'); title('Acc')
% figure; plot(gyrRaw.utcTimeMillis, [gyrRaw.UncalGyroXRadPerSec gyrRaw.UncalGyroYRadPerSec gyrRaw.UncalGyroZRadPerSec], '.'); 
% legend('X', 'Y', 'Z'); title('Gyr')
% xlabel('Time (ms)'); ylabel('Ang. vel. (rad/s)');
% figure; plot(magRaw.utcTimeMillis, [magRaw.UncalMagXMicroT magRaw.UncalMagYMicroT magRaw.UncalMagZMicroT], '.'); 
% legend('X', 'Y', 'Z'); title('Mag')
% xlabel('Time (ms)'); ylabel('Mag. field (uT)');
% 
% figure; plot(accRaw.utcTimeMillis-accRaw.utcTimeMillis(1), ~isnan(accRaw.UncalAccelXMps2), '.', 'MarkerSize', 5)
% hold on
% plot(gyrRaw.utcTimeMillis-accRaw.utcTimeMillis(1), 1*(~isnan(gyrRaw.UncalGyroXRadPerSec)), '.', 'MarkerSize', 5)
% plot(magRaw.utcTimeMillis-accRaw.utcTimeMillis(1), 1*(~isnan(magRaw.UncalMagXMicroT)), '.', 'MarkerSize', 5)
% legend('Acc', 'Gyr', 'Mag');
% 
% dtAcc = diff(accRaw.utcTimeMillis);
% [~, i] = max(dtAcc); dtAcc(i) = [];
% dtGyr = diff(gyrRaw.utcTimeMillis);
% [~, i] = max(dtGyr); dtGyr(i) = [];
% dtMag = diff(magRaw.utcTimeMillis);
% [~, i] = max(dtMag); dtMag(i) = [];
% 
% figure; histogram(dtAcc); xlabel('Acc dt (ms)');
% figure; histogram(dtGyr); xlabel('Gyr dt (ms)');
% figure; histogram(dtGyr); xlabel('Mag dt (ms)');
% 
% figure; plot((accRaw.elapsedRealtimeNanos-gyrRaw.elapsedRealtimeNanos(1))/1e6, ~isnan(accRaw.UncalAccelXMps2), '|', 'MarkerSize', 50,'LineWidth',2)
% hold on
% plot((gyrRaw.elapsedRealtimeNanos-gyrRaw.elapsedRealtimeNanos(1))/1e6, 1*(~isnan(gyrRaw.UncalGyroXRadPerSec)), '|', 'MarkerSize', 50,'LineWidth',2)
% plot((magRaw.elapsedRealtimeNanos-gyrRaw.elapsedRealtimeNanos(1))/1e6, 1*(~isnan(magRaw.UncalMagXMicroT)), '|', 'MarkerSize', 50,'LineWidth',2)
% legend('Acc', 'Gyr', 'Mag');
% xlabel('Time from start (ms)')
% grid on;
% 
% dtAcc = diff(accRaw.elapsedRealtimeNanos);
% [~, i] = max(dtAcc); dtAcc(i) = [];
% dtGyr = diff(gyrRaw.elapsedRealtimeNanos);
% [~, i] = max(dtGyr); dtGyr(i) = [];
% dtMag = diff(magRaw.elapsedRealtimeNanos);
% [~, i] = max(dtMag); dtMag(i) = [];
% 
% figure; histogram(dtAcc); xlabel('Acc chipset dt (ns)');
% figure; histogram(dtGyr); xlabel('Gyr chipset dt (ns)');
% figure; histogram(dtGyr); xlabel('Mag chipset dt (ns)');



%% Filter
% return
if ~filter
    isInvalid = obs(:, GnssLogUtils.COL_C1) > 50e6 | ...
                obs(:, GnssLogUtils.COL_C1) < -5e6 | ...
                obs(:, GnssLogUtils.COL_C2) > 50e6 | ...
                obs(:, GnssLogUtils.COL_C2) < -5e6;
    fprintf('%.2f%% of the measurements have been filtered due to invalid pseudoranges\n', sum(isInvalid)/length(isInvalid)*100);
    obs(isInvalid, :) = [];
end

isGpsIdx = obs(:, 3) == GnssLogUtils.OBS_ID_GPS;
isGloIdx = obs(:, 3) == GnssLogUtils.OBS_ID_GLO;
isGalIdx = obs(:, 3) == GnssLogUtils.OBS_ID_GAL;
isBdsIdx = obs(:, 3) == GnssLogUtils.OBS_ID_BDS;

%% Plots from RINEX
% indRnxGps = obs_rinex(:, 3) == GnssLogUtils.OBS_ID_GPS;
% indRnxBds = obs_rinex(:, 3) == GnssLogUtils.OBS_ID_BDS;
% indRnxGal = obs_rinex(:, 3) == GnssLogUtils.OBS_ID_GAL;

% constCount = [];
% for iRow = 1:size(obs_rinex)
%     switch obs_rinex(iRow, 3)
%         case 1 % GPS
%             if ~isnan(obs_rinex(iRow, 5))
%                 constCount(end+1) = GnssLogUtils.ID_GPS;
%             end
%             if ~isnan(obs_rinex(iRow, 9))
%                 constCount(end+1) = GnssLogUtils.ID_GPS;
%             end
%         case 3 % BDS
%             if ~isnan(obs_rinex(iRow, 5))
%                 constCount(end+1) = GnssLogUtils.ID_BDS;
%             end
%         case 4 % GAL
%             if ~isnan(obs_rinex(iRow, 5))
%                 constCount(end+1) = GnssLogUtils.ID_GAL;
%             end
%             if ~isnan(obs_rinex(iRow, 9))
%                 constCount(end+1) = GnssLogUtils.ID_GAL;
%             end
%     end
% end

% figure; histogram(constCount);
% xticklabels(GnssLogUtils.CONSTELLATIONS);
% xlabel('Constellation'); ylabel('# obs');
% saveas(gcf, [figsPath 'rnx_histogram_const.png'])
% 
% plotObsRnx(obs_rinex, obs_rinex, 5, isGpsIdx, 'RINEX GPS L1', 'pr (m)')
% plotObsRnx(obs_rinex, obs_rinex, 9, isGpsIdx, 'RINEX GPS L5', 'pr (m)')
% 
% plotObsRnx(obs_rinex, obs_rinex, 5, isGalIdx, 'RINEX Gal L1', 'pr (m)')
% plotObsRnx(obs_rinex, obs_rinex, 9, isGalIdx, 'RINEX Gal L5', 'pr (m)')

% close all;

%% Comparison with RINEX
nObsRnx = size(obs_rinex, 1);
% difPrGps = []; noObsGps = []; nGps = 0;
% difPrGal = []; noObsGal = []; nGal = 0;
% difPrBds = []; noObsBds = []; nBds = 0;

obsRnxConst = obs_rinex(indRnxGal, :); % indRnxGps
% svns = unique(obsRnxConst(:, 4));

obsLogConst = obs(isGalIdx, :); % isGpsIdx
svns = unique(obsLogConst(:, 4));

indDerConst = Derived.MEAS.ConstellationType == GnssLogUtils.ID_GAL;
obsDerConst = Derived.MEAS(indDerConst, :);

selFreq = GnssLogUtils.GPS_FREQUENCIES(2);

obsName = 'C'; % C L D S

switch obsName
    case 'C', colRnx = GnssLogUtils.COL_C2; units = '(m)';
        obsDerAll = obsDerConst.SmPrM;
    case 'L', colRnx = GnssLogUtils.COL_L2; units = '(cyc)';
        obsDerAll = obsDerConst.AdrM .* obsDerConst.CarrierFrequencyHz / LIGHTSPEED;
    case 'D', colRnx = GnssLogUtils.COL_D2; units = '(Hz)';
        obsDerAll = -obsDerConst.PrrMps .* obsDerConst.CarrierFrequencyHz / LIGHTSPEED;
    case 'S', colRnx = GnssLogUtils.COL_S2; units = '(dBHz)';
        obsDerAll = obsDerConst.Cn0DbHz;
end

% svns = [10 21];
for iSvn = 1:length(svns)
    % RINEX
    indRnxSvn = obsRnxConst(:, GnssLogUtils.COL_SVN) == svns(iSvn);
    towRnx = obsRnxConst(indRnxSvn, GnssLogUtils.COL_TOW);
    obsRnx = obsRnxConst(indRnxSvn, colRnx);
    
    % Processed GnssLog
    indObsSvn = obsLogConst(:, GnssLogUtils.COL_SVN) == svns(iSvn);
    towLog = obsLogConst(indObsSvn, GnssLogUtils.COL_TOW);
    obsLog = obsLogConst(indObsSvn, colRnx);
     
    % Derived data from GnssAnalysis
    indDerSvn = obsDerConst.Svid == svns(iSvn) & ...
                obsDerConst.CarrierFrequencyHz == selFreq;
    timeSecDer = obsDerConst.TimeNanos(indDerSvn) * 1e-9;
    obsDer = obsDerAll(indDerSvn);
    
    if ~all(isnan(obsRnx)) && ~isempty(obsLog)
        % Obs comparison
        figure;
        hold on;
        plot(towRnx, obsRnx, '.');
        plot(timeSecDer, obsDer, '.');
        plot(towLog, obsLog, '.');
        xlabel('TOW (s)'); ylabel([obsName units]);
        legend('RINEX', 'GnssAnalysis', 'GnssLog');
%         legend('GnssAnalysis', 'GnssLog');
        title(sprintf('SVN: %s', string(svns(iSvn))))
        hold off;
        
        % Obs difference
        indMeasRnx = ~isnan(obsRnx);
        obsLogInterpRnx = nan(size(obsRnx));
        obsLogInterpRnx(indMeasRnx) = interp1(towLog, obsLog, towRnx(indMeasRnx));
        diffRnxLog = obsRnx - obsLogInterpRnx;
        
        indMeasDer = ~isnan(obsDer);
        obsLogInterpDer = nan(size(obsDer));
        obsLogInterpDer(indMeasDer) = interp1(towLog, obsLog, timeSecDer(indMeasDer));
        diffDerLog = obsDer - obsLogInterpDer;
        
%         if obsName == 'C'
%             % Pr difference
%             figure;
%             plot(towRnx, -diffRnxLog);
%             xlabel('TOW (s)'); ylabel('Pr_{Log}-Pr_{RNX}');
%             title(sprintf('SVN: %s', string(svns(iSvn))))
%         end
        
%         figure;
%         subplot(2,1,1);
% %         plot(towRnx, diffRnxLog, '.');                    % Time series
% %         xlabel('TOW (s)'); ylabel([obsName ' ' units]);
%         histogram(diffRnxLog)                               % Histogram
%         xlabel([obsName ' ' units]);
%         title(sprintf('Error wrt RINEX SVN: %s', string(svns(iSvn))))
%         
%         subplot(2,1,2);
% %         plot(timeSecDer, diffDerLog, '.');                    % Time series
% %         xlabel('TOW (s)'); ylabel([obsName ' ' units]);
%         histogram(diffDerLog)                               % Histogram
%         xlabel([obsName ' ' units]);
%         title(sprintf('Error wrt GnssAnalysis SVN: %s', string(svns(iSvn))))
        
%         % CMC
%         cmcRinex = obsRnxConst(indRnxSvn, GnssLogUtils.COL_C1) - ...
%             obsRnxConst(indRnxSvn, GnssLogUtils.COL_L1) * LIGHTSPEED / GnssLogUtils.GPS_FREQUENCIES(1);
%         cmcGnssLog = obsLogConst(indObsSvn, GnssLogUtils.COL_C1) - ...
%             obsLogConst(indObsSvn, GnssLogUtils.COL_L1) * LIGHTSPEED / GnssLogUtils.GPS_FREQUENCIES(1);
%         
%         cmcRinex = cmcRinex - mean(cmcRinex, 'omitnan'); % Remove ambiguity
%         cmcGnssLog = cmcGnssLog - mean(cmcGnssLog, 'omitnan'); % Remove ambiguity
%         
%         % Remove very large uncertainties
%         indLarge = obsLogConst(:, GnssLogUtils.COL_U1L) > 1e30;
%         obsLogConst(indLarge, GnssLogUtils.COL_U1L) = 1;
%         
%         figure; hold on;
% %         plot(towRnx, cmcRinex, 'b.');
%         subplot(2,1,1)
%         plot(towLog, cmcGnssLog, 'r.');
%         xlabel('TOW'); ylabel('CMC (m)');
%         %         legend('GnssLog');
% %         legend('RINEX', 'GnssLog');
%         title(sprintf('SVN: %s', string(svns(iSvn))))
% %         subplot(2,1,2)
% % %         plot(towLog, 10*log10(obsLogConst(indObsSvn, GnssLogUtils.COL_U1L)), '.')
% %         semilogy(towLog, obsLogConst(indObsSvn, GnssLogUtils.COL_U1L), '.')
% %         grid on
% %         xlabel('TOW'); ylabel('\sigma_L (cyc)');
%         subplot(2,1,2)
%         plot(towLog, obsLogConst(indObsSvn, GnssLogUtils.COL_L1L), '.')
%         ylim([-0.5 1.5])
%         grid on
%         xlabel('TOW'); ylabel('LLI');
%         saveas(gcf, ['figs/' campaignName '/' phoneName '/cmc/lli_' num2str(svns(iSvn))], 'png');
    end
end

%% Functions
function plotMeas(obs, obsPlot, obsIdx, isConstIdx, freq, ttl, ylab)
    global figsPath 
    isFreqIdx = isConstIdx & (obs(:, 5) == freq);
    prns = unique(obs(isFreqIdx, 4));
    leg = cell(length(prns), 1);
    figure; hold on
    for iPrn = 1:length(prns)
        idx = isFreqIdx & (obs(:, 4) == prns(iPrn));
        plot(obs(idx, 2), obsPlot(idx, obsIdx), '.');
        leg{iPrn} = sprintf('prn: %d', prns(iPrn, 1));
    end
    legend(leg, 'Location', 'northeastoutside');
    xlabel('TOW (s)'); ylabel(ylab)
    title(ttl)
    saveas(gcf, [figsPath ttl '-' ylab '.png'])
end

function plotObsRnx(obs, obsPlot, obsIdx, isConstIdx, ttl, ylab)
    global figsPath 
    prns = unique(obs(isConstIdx, 4));
    leg = cell(length(prns), 1);
    figure; hold on
    for iPrn = 1:length(prns)
        idx = isConstIdx & (obs(:, 4) == prns(iPrn));
        plot(obs(idx, 2), obsPlot(idx, obsIdx), '.');
        leg{iPrn} = sprintf('prn: %d', prns(iPrn, 1));
    end
    legend(leg, 'Location', 'northeastoutside');
    xlabel('TOW (s)'); ylabel(ylab)
    title(ttl)
    saveas(gcf, [figsPath ttl '-' ylab '.png'])
end
