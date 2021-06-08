classdef GnssLogUtils < handle
    % GNSSUTILS Class with constants and utilities for processing GnssLog
    %   For more details about GnssLog flags and fields see
    %   https://developer.android.com/reference/android/location/GnssMeasurement.html
    %   and ESA's white paper "USING GNSS RAW MEASUREMENTS ON ANDROID DEVICES"
    %
    %     Author: Arnau Ochoa Banuelos (CS Group), March 2021
    
    %% Properties
    properties (Constant)
        CELERITY = 299792458;
        
        % Frequency bands
        GPS_BANDNAMES = {'L1' 'L5'};
        GPS_FREQUENCIES = [1575.42e6 1176.45e6];
        GLO_BANDNAMES = {'L1'};
        BDS_BANDNAMES = {'C2I'}
        BDS_FREQUENCIES = [1561.098e6];
        GAL_BANDNAMES = {'E1' 'E5a'};
        GAL_FREQUENCIES = [1575.42e6 1176.45e6];
        
        % Time relationships
        NUMBER_NANOSECONDS_HOUR     = 3600e9;
        NUMBER_NANOSECONDS_DAY      = 86400e9;
        NUMBER_NANOSECONDS_WEEK     = 604800e9;
        NUMBER_NANOSECONDS_100MILLI = 1e8;
        
        % Quality checks
%         MAX_BIAS_SECONDS = 10;      % Max bias for week rollover correction
        
        % Constellation id's (GnssLog)
        CONSTELLATIONS = {'G' 'S' 'R' 'Z' 'C' 'E' '-' '-' 'U'};
        ID_GPS = 1;
        ID_SBS = 2;
        ID_GLO = 3;
        ID_QZS = 4;
        ID_BDS = 5;
        ID_GAL = 6;
        ID_UNK = 9;
        
        % Constellation id's (obs_rinex)
        OBS_CONSTELLATIONS = {'G' 'R' 'C' 'E'};
        OBS_ID_GPS = 1; 
        OBS_ID_GLO = 2; 
        OBS_ID_BDS = 3; 
        OBS_ID_GAL = 4; 
        OBS_ID_SBS = 100;
        
        % Columns for obs_rinex matrix
        COL_WN = 1; COL_TOW = 2; COL_CONST = 3; COL_SVN = 4;        % Sat obs info
        COL_C1 = 5; COL_L1 = 6; COL_D1 = 7; COL_S1 = 8;             % 1st band obs
        COL_C2 = 9; COL_L2 = 10; COL_D2 = 11; COL_S2 = 12;          % 2nd band obs
        COL_U1C = 13; COL_U1L = 14; COL_U1D = 15; COL_L1L = 16;     % 1st band extra
        COL_U2C = 17; COL_U2L = 18; COL_U2D = 19; COL_L2L = 20;     % 2nd band extra
        COL_MAX = 20;
        
        % Observation types: 1 obs -> 3 characters
        % UfC: Code pseudorange uncertainty for frequency band f
        % UfL: Carrier phase uncertainty for frequency band f
        % UfD: Doppler uncertainty for frequency band f
        % LfL: Loss of Lock Indicator (LLI) for frequency band f
        OBS_TYPES_GPS = 'C1CL1CD1CS1CC5XL5XD5XS5XU1CU1LU1DL1LU5CU5LU5DL5L';
        OBS_TYPES_GLO = 'C1CL1CD1CS1C------------U1CU1LU1DL1L------------';
        OBS_TYPES_BDS = 'C2XL2XD2XS2X------------U2CU2LU2DL2L------------';
        OBS_TYPES_GAL = 'C1XL1XD1XS1XC5XL5XD5XS5XU1CU1LU1DL1LU5CU5LU5DL5L';
       
        % Loss of Lock Indicator Bit
        LLI_OK = 0;
        LLI_LOSSLOCK = 1;
        LLI_HALFCYCLE = 2;
    end
    
    %% Public methods
    methods (Static)
        function gnssRaw = filterFields(gnssRaw, validMeas)
            fn = fieldnames(gnssRaw);
            for iField = 1:length(fn)
                gnssRaw.(fn{iField}) = gnssRaw.(fn{iField})(validMeas);
            end
        end 
        
        function [trueFreq] = getTrueFrequencies(measFreq, refFreq)
            % Find indices of refFreq for the true frequencies closest to the measured frequencies
            [~, freqInd] = min(abs(measFreq - refFreq), [], 2);
            % Correct frequencies
            trueFreq = refFreq(freqInd);
        end 
        
        function idObsConst = getIdObsConst(idConst)
            % Transform from GnssLog's constellation id to obs_rinex's constellation id
            switch idConst
                case GnssLogUtils.ID_GPS, idObsConst = GnssLogUtils.OBS_ID_GPS;
                case GnssLogUtils.ID_SBS, idObsConst = GnssLogUtils.OBS_ID_SBS;
                case GnssLogUtils.ID_GLO, idObsConst = GnssLogUtils.OBS_ID_GLO;
                case GnssLogUtils.ID_BDS, idObsConst = GnssLogUtils.OBS_ID_BDS;
                case GnssLogUtils.ID_GAL, idObsConst = GnssLogUtils.OBS_ID_GAL;
                otherwise, idObsConst = -1;
            end
        end
        
        function constIds = getIdsObsConstFromStr(constStr)
            % Transform from constellation letters to constellation id's
            % used in rinex obs matrix
            constIds = nan(size(constStr));
            for iConst = 1:length(constStr)
                constIds(iConst) = find(strcmp(GnssLogUtils.OBS_CONSTELLATIONS, constStr(iConst)));
            end
        end

        function svnIndSorted = sortIndices(svnInd, frequencies, idConst)
            % Sort indices to match with desired order (set in GnssLogUtils)
            switch idConst
                case GnssLogUtils.ID_GPS, [~, idx] = intersect(frequencies, GnssLogUtils.GPS_FREQUENCIES);
                case GnssLogUtils.ID_BDS, [~, idx] = intersect(frequencies, GnssLogUtils.BDS_FREQUENCIES);
                case GnssLogUtils.ID_GAL, [~, idx] = intersect(frequencies, GnssLogUtils.GAL_FREQUENCIES);
                otherwise, error('This constellation is not supported')
            end
            svnIndSorted = svnInd(idx);
        end
    end
end