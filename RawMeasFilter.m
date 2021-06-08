classdef RawMeasFilter < handle
    % RAWMEASFILTER Class for filtering Android raw GNSS measurements
    %   For more details about GnssLog flags and fields see
    %   https://developer.android.com/reference/android/location/GnssMeasurement.html
    %   and ESA's white paper "USING GNSS RAW MEASUREMENTS ON ANDROID DEVICES"
    %
    %     Author: Arnau Ochoa Banuelos (CS Group), March 2021
    
    %% Properties
    properties (Constant)
        % Time states
        STATE_TOW_KNOWN = 16384;                % If TOW decoded -> TOW known
        STATE_GLO_TOD_KNOWN = 32768;            % If TOD decoded -> TOD known
        STATE_GAL_E1C_2ND_CODE_LOCK = 2048;     % Recommended by ESA's white paper
        
        % ADR states
        ADR_STATE_VALID = 1;                    % State of ADR is valid
        ADR_STATE_RESET = 2;                    % Reset of ADR detected
        ADR_STATE_CYCLE_SLIP = 4;               % Cycle slip detected
        ADR_STATE_HALF_CYCLE_RESOLVED = 8;      % Reports whether the half cycle ambiguity has been resolved.
        ADR_STATE_HALF_CYCLE_REPORTED = 16;     % Reports whether a half cycle has been reported by the GNSS hardware.
        
        % Thresholds
        THRES_BIAS_UNCERTAINTY_NANOS = 1e8;     % Quality indicator of rx clock
        THRES_RX_SVT_UNCERTAINTY_NANOS = 1e8;   % Received sv tx time uncertainty
        THRES_PRR_UNCERTAINTY_MPSEC = 1e8;      % Doppler uncertainty
        THRES_ADR_UNCERTAINTY_CYC = 1e8;        % Carrier phase uncertainty
        
    end
    
    %% Public methods
    methods (Static)
        function valid = checkTimeState(states, constTypes)
            % Returns boolean array of valid states considering time
            % related flags
            
            valid = zeros(size(states));
            % GPS
            gpsStates = zeros(size(states));
            gpsStates(constTypes == GnssLogUtils.ID_GPS) = states(constTypes == GnssLogUtils.ID_GPS);
            valid = valid | RawMeasFilter.checkGpsTime(gpsStates);
            % GLO
            gloStates = zeros(size(states));
            gloStates(constTypes == GnssLogUtils.ID_GLO) = states(constTypes == GnssLogUtils.ID_GLO);
            valid = valid | RawMeasFilter.checkGloTime(gloStates);
            % GAL
            galStates = zeros(size(states));
            galStates(constTypes == GnssLogUtils.ID_GAL) = states(constTypes == GnssLogUtils.ID_GAL);
            valid = valid | RawMeasFilter.checkGalTime(galStates);
            % BDS
            bdsStates = zeros(size(states));
            bdsStates(constTypes == GnssLogUtils.ID_BDS) = states(constTypes == GnssLogUtils.ID_BDS);
            valid = valid | RawMeasFilter.checkBdsTime(bdsStates);
        end
        
        function valid = findConstellations(constTypes, constToKeep)
            valid = zeros(size(constTypes));
            for iConst = 1:length(constToKeep)
                valid = valid | constTypes == constToKeep(iConst);
            end
        end
        
        function valid = checkAdrState(states)
            % Returns boolean array of valid states considering ADR
            % (i.e. carrier phase) related flags
            valid = RawMeasFilter.haveState(states, RawMeasFilter.ADR_STATE_VALID) ...  % ADR_VALID == 1
                & ~RawMeasFilter.haveState(states, RawMeasFilter.ADR_STATE_RESET) ...   % ADR_RESET == 0
                & ~RawMeasFilter.haveState(states, RawMeasFilter.ADR_STATE_CYCLE_SLIP); % ADR_CYCLE_SLIP == 0
        end
        
        function valid = checkBounded(gnssRaw)
            % Returns boolean array of valid measurements considering
            % measurements that are bounded or can have invalid values
            
            valid = ones(size(gnssRaw.State));
            % Constellation can't be unknown
            valid = valid & gnssRaw.ConstellationType ~= GnssLogUtils.ID_UNK;
            % TimeNanos can't be empty
            valid = valid & ~isnan(gnssRaw.TimeNanos);
            % FullBiasNanos can't be 0 or empty
            valid = valid & (gnssRaw.FullBiasNanos ~= 0 & ~isnan(gnssRaw.FullBiasNanos));
            % BiasUncertaintyNanos can't be too large
            valid = valid & (gnssRaw.BiasUncertaintyNanos ...
                < RawMeasFilter.THRES_BIAS_UNCERTAINTY_NANOS);
            % ReceivedSvTimeUncertaintyNanos can't be too large
            valid = valid & (gnssRaw.ReceivedSvTimeUncertaintyNanos ...
                < RawMeasFilter.THRES_RX_SVT_UNCERTAINTY_NANOS);
            % PseudorangeRateUncertaintyMetersPerSecond (Doppler) can't be too large
%             valid = valid & (gnssRaw.PseudorangeRateUncertaintyMetersPerSecond ...
%                 < RawMeasFilter.THRES_PRR_UNCERTAINTY_MPSEC);
            % AccumulatedDeltaRangeUncertaintyMeters (Carrier phase) can't be too large
%             valid = valid & (gnssRaw.AccumulatedDeltaRangeUncertaintyMeters ./ ...
%                     (GnssLogUtils.CELERITY./gnssRaw.CarrierFrequencyHz) ...
%                     < RawMeasFilter.THRES_ADR_UNCERTAINTY_CYC);
        end
    end
    
    %% Private methods
    methods (Static, Access = private)
        function valid = checkGpsTime(states)
            % Check time state for GPS
            valid = RawMeasFilter.haveState(states, RawMeasFilter.STATE_TOW_KNOWN);
        end
        function valid = checkGloTime(states)
            % Check time state for GLONASS
            valid = RawMeasFilter.haveState(states, RawMeasFilter.STATE_GLO_TOD_KNOWN);
        end
        function valid = checkGalTime(states)
            % Check time state for GALILEO
            valid = RawMeasFilter.haveState(states, RawMeasFilter.STATE_TOW_KNOWN); %& ...
            %    RawMeasFilter.haveState(states, RawMeasFilter.STATE_GAL_E1C_2ND_CODE_LOCK);
        end
        function valid = checkBdsTime(states)
            % Check time state for BEIDOU
            valid = RawMeasFilter.haveState(states, RawMeasFilter.STATE_TOW_KNOWN);
        end
        function indices = haveState(states, required)
            indices = bitand(states, required) == required;
        end
    end
end