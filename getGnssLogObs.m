function [obsRinex, obsRinexType, obsRinexUtcMillis, gnssAnalysis, accMeas, gyrMeas, magMeas, sensorAnalysis] ...
    = getGnssLogObs(dirname, filename, filter)
% GETGNSSLOGOBS Reads the provided GnssLog and process the data to output
% the GNSS observables in RINEX format as well as the INS+Mag measurements
%
%   [obs, obsType, gnssAnalysis, accMeas, gyrMeas, magMeas, sensorAnalysis] ...
%       = GETGNSSLOGOBS(dirname, filename, filter)
%
% Input:
%   dirName             = string with directory of fileName
%   fileName            = string with filename
%   filter              = boolean flag to activate/deactivate filtering of measurements.
%
% Output:  
%   obsRinex            = matrix containing observations and extra info form GnssLog
%   obsRinexType        = description of the observation types by constellation
%   obsRinexUtcMillis   = vector of timestamps of obsRinex epochs in UTC
%   gnssAnalysis        = structure containing an analysis of gnssRaw (e.g. missing fields)
%   accMeas             = structure containing accelerometer measurements from log file
%   gyrMeas             = structure containing gyroscope measurements from log file
%   magMeas             = structure containing magnetometer measurements from log file
%   sensorAnalysis      = structure containing an analysis of accMeas, gyrMeas and
%                           magMeas (e.g. missing fields)
%
% Dependencies: This function uses DataHash, see: 
%       https://fr.mathworks.com/matlabcentral/fileexchange/31272-datahash

% Author: Arnau Ochoa Banuelos (CS Group), April 2021

inputSignature = DataHash([dirname filename], 'file');
if Fcache.hasValidCacheFile(mfilename, inputSignature)
    disp('Loading Android raw measurements...');
    cache = Fcache.get(mfilename, inputSignature);
    obsRinex = cache.obsRinex;
    obsRinexType = cache.obsRinexType;
    obsRinexUtcMillis = cache.obsRinexUtcMillis;
    gnssAnalysis = cache.gnssAnalysis;
    accMeas = cache.accMeas;
    gyrMeas = cache.gyrMeas;
    magMeas = cache.magMeas;
    sensorAnalysis = cache.sensorAnalysis;
    return;
else
    disp([mfilename '>> No cached file found, reading GnssLog file']);
end


disp('Reading Android raw measurements...');
[gnssMeas, gnssAnalysis, accMeas, gyrMeas, magMeas, sensorAnalysis] = readGnssLog(dirname, filename);
disp('Processing Android GNSS raw measurements...');
[obsRinex, obsRinexType, obsRinexUtcMillis] = processGnssRaw(gnssMeas, filter);
% Check obsRinexUtcMillis
[uniqueUtc, iUtcA, ~] = unique(obsRinexUtcMillis);
[uniqueObs, iObsA, ~] = unique(obsRinex(:, 2));
assert(length(uniqueUtc) == length(uniqueObs),'obsRinexUtcTimeMillis does not have the same number of epochs as obsRinex');
assert(all(iObsA == iUtcA), 'Correspondecies in obsRinexUtcTimeMillis do not match tow in obsRinex');

cache = struct();
cache.obsRinex = obsRinex;
cache.obsRinexType = obsRinexType;
cache.obsRinexUtcMillis = obsRinexUtcMillis;
cache.gnssAnalysis = gnssAnalysis;
cache.accMeas = accMeas;
cache.gyrMeas = gyrMeas;
cache.magMeas = magMeas;
cache.sensorAnalysis = sensorAnalysis;
Fcache.cache(mfilename, inputSignature, cache);
end