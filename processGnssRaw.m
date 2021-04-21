function [obs, obsType] = processGnssRaw(gnssRaw, filter)
% PROCESSGNSSRAW computes the RINEX observations from Android GNSS raw measurements
%
%   [obs, obsType] = PROCESSGNSSRAW(gnssRaw, filter)
% 
% Input:  
%   gnssRaw = structure containing the Android GNSS raw measurements (from readGnssLog).
%   filter = boolean flag to activate/deactivate filtering of measurements.
% Output:  
%   obs = matrix containing observations and extra info form GnssLog
%   obsType = description of the observation types by constellation

% Author: Arnau Ochoa Banuelos (CS Group), March 2021

% factored into a few main sub-functions:
% correctFrequencies()
% computeSystemTimes()
% checkWeekRollover()
% filterRawMeas()
% packObservations()

%% Keep only desired constellations
desiredConst = [GnssLogUtils.ID_GPS GnssLogUtils.ID_GAL GnssLogUtils.ID_BDS];
isDesiredInd = RawMeasFilter.findConstellations(gnssRaw.ConstellationType, desiredConst);
gnssRaw = GnssLogUtils.filterFields(gnssRaw, isDesiredInd);
clear isDesiredInd desiredConst

%% Correct frequencies from measured to true
gnssRaw = correctFrequencies(gnssRaw);

%% Filter raw measurements
gnssRaw = filterRawMeas(gnssRaw, filter);

%% Obtain reception time in system time for each constellation
% Compute arrival time wrt reference time (GPS) using gnssRaw.FullBiasNanos(1), 
% so that tRxNanos includes rx clock drift since the first epoch:
tRxRefNanos = gnssRaw.TimeNanos - gnssRaw.FullBiasNanos(1); % Remove first FullBias bc both are int64

% Obtain reception times relative to each system time (e.g. GPST, GLONASST)
tRxSysNanos = computeSystemTimes(tRxRefNanos, gnssRaw);

% Subtract fractional part of bias and offset after modulo operation
tRxSysNanos = tRxSysNanos - gnssRaw.TimeOffsetNanos - gnssRaw.BiasNanos;
tRxSysSec = tRxSysNanos*1e-9;
% clear tRxSysNanos

%% Filter measurements by arrival time
% Reception time can't be negative
validMeas = tRxSysSec >= 0;
tRxRefNanos = tRxRefNanos(validMeas);
tRxSysSec = tRxSysSec(validMeas);
gnssRaw = GnssLogUtils.filterFields(gnssRaw, validMeas);
fprintf('\t%.2f%% of the measurements have been filtered due to negative TOA\n', ...
        100*(length(validMeas)-sum(validMeas))/length(validMeas));
clear validMeas

%% Compute code pseudoranges
% Transmission time wrt each system time
tTxSec = double(gnssRaw.ReceivedSvTimeNanos)*1e-9;
% Pseudorange in seconds
prSeconds = (tRxSysSec - tTxSec);
% Correct rollover
[prSeconds, tRxRefNanos] = checkWeekRollover(prSeconds, tRxRefNanos);
% Pseudorange in meters
prMeas = prSeconds.* MagnitudeConstants.celerity;

%% Compute carrier phase measurements
% Tranform from meters to cycle units
carrierMeasCyc = gnssRaw.AccumulatedDeltaRangeMeters .* ...
                gnssRaw.CarrierFrequencyHz ./ MagnitudeConstants.celerity;

%% Compute Doppler measurements
doppMeasHz = -gnssRaw.PseudorangeRateMetersPerSecond .* ...
                gnssRaw.CarrierFrequencyHz ./ MagnitudeConstants.celerity;

%% Compute pseudorange uncertainty
prSigmaNanos = sqrt(gnssRaw.TimeUncertaintyNanos.^2 + ...
                    gnssRaw.BiasUncertaintyNanos.^2 + ...
                    double(gnssRaw.ReceivedSvTimeUncertaintyNanos).^2);
prSigmaMeters = prSigmaNanos .* 1e-9 .* MagnitudeConstants.celerity;

%% Compute carrier phase uncertainty
carrierSigmaCyc = gnssRaw.AccumulatedDeltaRangeUncertaintyMeters .* ...
                gnssRaw.CarrierFrequencyHz ./ MagnitudeConstants.celerity;

%% Compute Doppler uncertainty
dopplerSigmaHz = gnssRaw.PseudorangeRateUncertaintyMetersPerSecond .* ...
                gnssRaw.CarrierFrequencyHz ./ MagnitudeConstants.celerity;

%% Compute loss of lock indicator bit
lliBit = nan(size(carrierMeasCyc));
% Carrier meas is ok
isAdrOk = RawMeasFilter.checkAdrState(gnssRaw.AccumulatedDeltaRangeState);
lliBit(isAdrOk) = GnssLogUtils.LLI_OK;
% Loss of Lock/Cycle Slip
isAdrCS = gnssRaw.AccumulatedDeltaRangeState == RawMeasFilter.ADR_STATE_CYCLE_SLIP;
lliBit(isAdrCS) = GnssLogUtils.LLI_LOSSLOCK;
% Half cycle slip/Opposite wavelength factor
isAdrHCS = gnssRaw.AccumulatedDeltaRangeState == RawMeasFilter.ADR_STATE_HALF_CYCLE_REPORTED;
% nansum to account for possible combination of CS and HCS, see LLI specs
lliBit(isAdrHCS) = nansum([lliBit(isAdrHCS) GnssLogUtils.LLI_HALFCYCLE*ones(size(lliBit(isAdrHCS)))], 2); 

%% Compute signal strength indicator bit (?)


%% Arrange measurements in obs matrix
wNum = floor(-double(gnssRaw.FullBiasNanos)/GnssLogUtils.NUMBER_NANOSECONDS_WEEK);
tow = double(mod(tRxRefNanos, GnssLogUtils.NUMBER_NANOSECONDS_WEEK))*1e-9;

[obs, obsType] = packObservations(wNum, tow, gnssRaw, prMeas, carrierMeasCyc, ...
                    doppMeasHz, prSigmaMeters, carrierSigmaCyc, dopplerSigmaHz, lliBit);

end %end of function processGnssRaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gnssRaw = correctFrequencies(gnssRaw)
    % GPS
    gpsMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GPS);
    gnssRaw.CarrierFrequencyHz(gpsMeasIdx) = ...
        GnssLogUtils.getTrueFrequencies(gnssRaw.CarrierFrequencyHz(gpsMeasIdx), GnssLogUtils.GPS_FREQUENCIES);
    % BDS
    bdsMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_BDS);
    gnssRaw.CarrierFrequencyHz(bdsMeasIdx) = ...
        GnssLogUtils.getTrueFrequencies(gnssRaw.CarrierFrequencyHz(bdsMeasIdx), GnssLogUtils.BDS_FREQUENCIES);
    % GAL
    galMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GAL);
    gnssRaw.CarrierFrequencyHz(galMeasIdx) = ...
        GnssLogUtils.getTrueFrequencies(gnssRaw.CarrierFrequencyHz(galMeasIdx), GnssLogUtils.GAL_FREQUENCIES);
    
end %end of function correctFrequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tRxSysNanos] = computeSystemTimes(tRxRefNanos, gnssRaw)
    % Perform modulo operation to obtain system time for each constellation
    tRxSysNanos = nan(size(tRxRefNanos));
    % GPS
    gpsMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GPS);
    tRxSysNanos(gpsMeasIdx) = double(mod(tRxRefNanos(gpsMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_WEEK));
    % GLONASS
    gloMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GLO);
    tRxRefNanos(gloMeasIdx) = mod(tRxRefNanos(gloMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_WEEK);
    tRxSysNanos(gloMeasIdx) = double(mod(tRxRefNanos(gloMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_DAY)) ...
                                - (3*GnssLogUtils.NUMBER_NANOSECONDS_HOUR) + gnssRaw.LeapSecond(gloMeasIdx);
    % GALILEO
    galMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GAL);
    tRxSysNanos(galMeasIdx) = double(mod(tRxRefNanos(galMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_WEEK)); 
    % GALILEO E1C 2nd code decoded
    gal2ndCodeMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_GAL & ...
                            gnssRaw.State == RawMeasFilter.STATE_GAL_E1C_2ND_CODE_LOCK); % check e1c 2nd code
    tRxSysNanos(gal2ndCodeMeasIdx) = double(mod(tRxRefNanos(gal2ndCodeMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_100MILLI));
    % BEIDOU
    bdsMeasIdx = find(gnssRaw.ConstellationType == GnssLogUtils.ID_BDS);
    tRxSysNanos(bdsMeasIdx) = double(mod(tRxRefNanos(bdsMeasIdx), GnssLogUtils.NUMBER_NANOSECONDS_WEEK)) - 14e9;
end %end of function computeSystemTimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prSeconds, tRxRefNanos] = checkWeekRollover(prSeconds, tRxRefNanos)
    isRollover = prSeconds > GnssLogUtils.NUMBER_NANOSECONDS_WEEK/2e9;
    if any(isRollover)
        fprintf('\t\nWeek rollover detected in time tags. Adjusting ...\n')
        prS = prSeconds(isRollover);
        delS = round(prS/(GnssLogUtils.NUMBER_NANOSECONDS_WEEK/1e9))*(GnssLogUtils.NUMBER_NANOSECONDS_WEEK/1e9);
        prS = prS - delS;
        % prS are in the range [-WEEKSEC/2 : WEEKSEC/2];

        % check that common bias is not huge
%         assert(abs(prS) < GnssLogUtils.MAX_BIAS_SECONDS, 'Failed to correct week rollover');
        prSeconds(isRollover) = prS; %put back into prSeconds vector
        % Now adjust tRxSeconds by the same amount:
        tRxRefNanos(isRollover) = tRxRefNanos(isRollover) - int64(delS*1e9);
        fprintf('\tCorrected week rollover\n')
    end
end %end of function checkWeekRollover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gnssRaw = filterRawMeas(gnssRaw, filter)
    if filter
        % Filter by state
        validMeasTimeState = RawMeasFilter.checkTimeState(gnssRaw.State, gnssRaw.ConstellationType);
        fprintf('\t%d out of %d will be removed due to invalid time state\n', ...
                length(validMeasTimeState)-sum(validMeasTimeState), length(validMeasTimeState));
%         validMeasAdrState = RawMeasFilter.checkAdrState(gnssRaw.AccumulatedDeltaRangeState);
        validMeasAdrState = ones(size(gnssRaw.TimeNanos));
        fprintf('\t%d out of %d will be removed due to invalid ADR state\n', ...
                length(validMeasAdrState)-sum(validMeasAdrState), length(validMeasAdrState));

        % Filter by values that have a threshold
        validMeasBounded = RawMeasFilter.checkBounded(gnssRaw);
        fprintf('\t%d out of %d will be removed due to invalid or bounded values\n', ...
                length(validMeasBounded)-sum(validMeasBounded), length(validMeasBounded));

        validMeas = validMeasTimeState & validMeasAdrState & validMeasBounded;

        % Keep only valid measurements on all fields of gnssRaw
        gnssRaw = GnssLogUtils.filterFields(gnssRaw, validMeas);
        fprintf('\t%.2f%% of the measurements have been filtered \n', ...
                100*(length(validMeas)-sum(validMeas))/length(validMeas));
        clear validMeas validMeasTimeState validMeasAdrState validMeasBounded
    end

    %check for fields that are commonly all zero and may be missing from gnssRaw
    if ~isfield(gnssRaw,'TimeUncertaintyNanos')
        gnssRaw.TimeUncertaintyNanos = zeros(size(gnssRaw.TimeNanos));
    end
    if ~isfield(gnssRaw,'BiasNanos')
        gnssRaw.BiasNanos = zeros(size(gnssRaw.TimeNanos));
    end
    if ~isfield(gnssRaw,'TimeOffsetNanos')
        gnssRaw.TimeOffsetNanos = zeros(size(gnssRaw.TimeNanos));
    end
    if ~isfield(gnssRaw,'LeapSecond')
        gnssRaw.LeapSecond = zeros(size(gnssRaw.TimeNanos));
    end
    % New fields that may not be available in some datasets
    if ~isfield(gnssRaw,'BasebandCn0DbHz')
        gnssRaw.BasebandCn0DbHz = nan(size(gnssRaw.TimeNanos));
    end
    if ~isfield(gnssRaw,'CodeType')
        gnssRaw.CodeType = nan(size(gnssRaw.TimeNanos));
    end
    if ~isfield(gnssRaw,'ChipsetElapsedRealtimeNanos')
        gnssRaw.ChipsetElapsedRealtimeNanos = nan(size(gnssRaw.TimeNanos));
    end
end %end of function filterRawMeas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obs, obsType] = packObservations(wNumVec, towVec, gnssRaw, prMeasVec, carrierMeasVec, ...
    doppMeasVec, prSigmaMetersVec, carrierSigmaCycVec, dopplerSigmaHzVec, lliBitVec)

    % Number of observations: unique combinations of wnum, tow, const and svid
    nObs = length(unique([wNumVec towVec gnssRaw.ConstellationType gnssRaw.Svid], 'stable', 'rows'));
    obs = nan(nObs, GnssLogUtils.COL_MAX);
    iObs = 0;
    
    wNumUnique = unique(wNumVec);
    for wNum = wNumUnique
        % Indices of current week num
        isWNumInd = wNumVec == wNum;
        % Unique tow in current wNum
        towUnique = unique(towVec(isWNumInd));
        for tow = towUnique'
            % Indices of current tow & wNum
            isTowInd = isWNumInd & towVec == tow;
            % Unique constellations in current tow
            constUnique = unique(gnssRaw.ConstellationType(isTowInd));
            for const = constUnique'
                % Indices of current const & TOW & wNum
                isConstInd = isTowInd & gnssRaw.ConstellationType == const;
                % Unique constellations in current tow & const
                svnUnique = unique(gnssRaw.Svid(isConstInd));
                for svn = svnUnique'
                    iObs = iObs + 1;
                    obs(iObs, GnssLogUtils.COL_WN) = wNum;
                    obs(iObs, GnssLogUtils.COL_TOW) = tow;
                    obs(iObs, GnssLogUtils.COL_CONST) = GnssLogUtils.getIdObsConst(const);
                    obs(iObs, GnssLogUtils.COL_SVN) = svn;
                    
                    % Indices of current svn & const & TOW & wNum
                    svnInd = find(isConstInd & gnssRaw.Svid == svn);
                    svnInd = GnssLogUtils.sortIndices(svnInd, gnssRaw.CarrierFrequencyHz(svnInd), const);      
                    for iMeas = svnInd'
                        % First frequency band (4 columns)
                        if  gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.GPS_FREQUENCIES(1) || ...
                            gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.BDS_FREQUENCIES(1) || ...
                            gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.GAL_FREQUENCIES(1)
                                obs(iObs, GnssLogUtils.COL_C1) = prMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_L1) = carrierMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_D1) = doppMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_S1) = gnssRaw.Cn0DbHz(iMeas);
                                obs(iObs, GnssLogUtils.COL_U1C) = prSigmaMetersVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_U1L) = carrierSigmaCycVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_U1D) = dopplerSigmaHzVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_L1L) = lliBitVec(iMeas);
                        % Second frequency band (4 columns)
                        elseif  gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.GPS_FREQUENCIES(2) || ...
                            gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.BDS_FREQUENCIES(2) || ...
                            gnssRaw.CarrierFrequencyHz(iMeas) == GnssLogUtils.GAL_FREQUENCIES(2)
                                obs(iObs, GnssLogUtils.COL_C2) = prMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_L2) = carrierMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_D2) = doppMeasVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_S2) = gnssRaw.Cn0DbHz(iMeas);
                                obs(iObs, GnssLogUtils.COL_U2C) = prSigmaMetersVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_U2L) = carrierSigmaCycVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_U2D) = dopplerSigmaHzVec(iMeas);
                                obs(iObs, GnssLogUtils.COL_L2L) = lliBitVec(iMeas);
                        end
                    end
                end
            end
        end
    end
    
    % Build obsType
    if mod(length(GnssLogUtils.OBS_TYPES_GPS), 3) ~= 0, error('OBS_TYPES_GPS is wrong'), end
    if mod(length(GnssLogUtils.OBS_TYPES_GLO), 3) ~= 0, error('OBS_TYPES_GLO is wrong'), end
    if mod(length(GnssLogUtils.OBS_TYPES_BDS), 3) ~= 0, error('OBS_TYPES_BDS is wrong'), end
    if mod(length(GnssLogUtils.OBS_TYPES_GAL), 3) ~= 0, error('OBS_TYPES_GAL is wrong'), end
    
    % GPS
    obsType(GnssLogUtils.OBS_ID_GPS).num_obs = length(GnssLogUtils.OBS_TYPES_GPS)/3;
    obsType(GnssLogUtils.OBS_ID_GPS).type_str = GnssLogUtils.OBS_TYPES_GPS;
    % GLONASS
    obsType(GnssLogUtils.OBS_ID_GLO).num_obs = length(GnssLogUtils.OBS_TYPES_GLO)/3;
    obsType(GnssLogUtils.OBS_ID_GLO).type_str = GnssLogUtils.OBS_TYPES_GLO;
    % BEIDOU
    obsType(GnssLogUtils.OBS_ID_BDS).num_obs = length(GnssLogUtils.OBS_TYPES_BDS)/3;
    obsType(GnssLogUtils.OBS_ID_BDS).type_str = GnssLogUtils.OBS_TYPES_BDS;
    % GALILEO
    obsType(GnssLogUtils.OBS_ID_GAL).num_obs = length(GnssLogUtils.OBS_TYPES_GAL)/3;
    obsType(GnssLogUtils.OBS_ID_GAL).type_str = GnssLogUtils.OBS_TYPES_GAL;

end %end of function packObservations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
