function [gnssMeas, gnssAnalysis, accMeas, gyrMeas, magMeas, sensorAnalysis] = readGnssLog(dirName,fileName)
% READGNSSLOG Read the log file created by Gnss Logger App in Android
%
%   [gnssMeas, gnssAnalysis, accMeas, gyrMeas, magMeas, sensorAnalysis] = READGNSSLOG(dirName,fileName)
%
% Input:
%   dirName         = string with directory of fileName,
%                   e.g. '~/Documents/MATLAB/Pseudoranges/2016-03-28'
%   fileName        = string with filename
%
% Output:
%   gnssMeas         = all GnssClock and GnssMeasurement fields from log file, including:
%                       .TimeNanos (int64)
%                       .FullBiasNanos (int64)
%                       ...
%                       .Svid
%                       .ReceivedSvTimeNanos (int64)
%                       .PseudorangeRateMetersPerSecond
%                       ...
%   gnssAnalysis    = structure containing an analysis of gnssRaw (e.g. missing fields)
%   accMeas         = structure containing accelerometer measurements from log file
%   gyrMeas         = structure containing gyroscope measurements from log file
%   magMeas         = structure containing magnetometer measurements from log file
%   sensorAnalysis  = structure containing an analysis of accMeas, gyrMeas and
%                       magMeas (e.g. missing fields)

% Author: Arnau Ochoa Banuelos (CS Group), March 2021
% Based on readGnssLog function by Frank van Diggelen. Modified so that new 
% GNSS data types can be read (e.g. chars in CodeType) as well as
% Accelerometer, Gyroscope and Magnetometer measurements.
% (see github.com/google/gps-measurement-tools to see the original function)

% factored into a few main sub-functions:
% MakeCsv()             *modified from readGnssLog*
% ReadRawCsv()          *modified from readGnssLog*
% PackGnssRaw()         *modified from readGnssLog*
% PackMeas()            *new*
% CheckGnssClock()
% ReportMissingFields() *modified from readGnssLog*

%% Initialize outputs and inputs
gnssAnalysis.GnssClockErrors = 'GnssClock Errors.';
gnssAnalysis.GnssMeasurementErrors = 'GnssMeasurement Errors.';
gnssAnalysis.ApiPassFail = '';

%% check we have the right kind of fileName
extension = fileName(end-3:end);
if ~any(strcmp(extension,{'.txt','.csv'}))
    error('Expecting file name of the form "*.txt", or "*.csv"');
end

%% read log file into a numeric matrix 'S', and a cell array 'header'
fullLogFilePath = [dirName,fileName];
gnssCsvFilePath = [dirName, fileName(1:end-4), '_Gnss.csv'];
accCsvFilePath = [dirName, fileName(1:end-4), '_UncalAccel.csv'];
gyrCsvFilePath = [dirName, fileName(1:end-4), '_UncalGyro.csv'];
magCsvFilePath = [dirName, fileName(1:end-4), '_UncalMag.csv'];
if ~isfile(gnssCsvFilePath) % Check if csv has already been created
    MakeCsv(fullLogFilePath, gnssCsvFilePath, accCsvFilePath, gyrCsvFilePath, magCsvFilePath);
end

[rawDataTable] = readRawCsv(gnssCsvFilePath); % GNSS raw
[accDataTable] = readRawCsv(accCsvFilePath); % Accelerometer
[gyrDataTable] = readRawCsv(gyrCsvFilePath); % Gyroscope
[magDataTable] = readRawCsv(magCsvFilePath); % Magnetometer

%% pack data into structure
[gnssMeas,missingGnss] = PackGnssRaw(rawDataTable);
[accMeas, accMissing] = PackMeas(accDataTable, 'acc');
[gyrMeas, gyrMissing] = PackMeas(gyrDataTable, 'gyr');
[magMeas, magMissing] = PackMeas(magDataTable, 'mag');

%% check clock and measurements
[gnssMeas,gnssAnalysis] = CheckGnssClock(gnssMeas,gnssAnalysis);
[gnssAnalysis, sensorAnalysis] = ReportMissingFields(gnssAnalysis,missingGnss, accMissing, gyrMissing, magMissing);

end %end of function ReadGnssLogger
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeCsv(fullLogFilePath, gnssCsvFilePath, accCsvFilePath, gyrCsvFilePath, magCsvFilePath)
%% make csv file, if necessary.
%And return extended csv file name (i.e. with full path in the name)

fprintf('\nReading file %s\n', fullLogFilePath)

%% read version
txtfileID = fopen(fullLogFilePath,'r');
if txtfileID<0
    error('file ''%s'' not found',fullLogFilePath);
end
line='';
while ~contains(lower(line),'version')
    line = fgetl(txtfileID);
    if ~ischar(line) %eof or error occurred
        if isempty(line)
            error('\nError occurred while reading file %s\n',fullLogFilePath)
        end
        break
    end
end
if line==-1
    fprintf('\nCould not find "Version" in input file %s\n',fullLogFilePath)
    return
end
%look for the beginning of the version number, e.g. 1.4.0.0
iDigits = regexp(line,'\d'); %index into the first number found in line
v = sscanf(line(iDigits(1):end),'%d.%d.%d.%d',4);
if length(v)<4
    v(end+1:4,1)=0; %make v into a length 4 column vector
end
%Now extract the platform
k = strfind(line,'Platform:');
if any(k)
    sPlatform = line(k+9:end);
else
    sPlatform = '';%set empty if 'Platform:' not found
end
if ~contains(sPlatform,'N') && ~contains(sPlatform,'10')
    %add || strfind(platform,'O') and so on for future platforms
    fprintf('\nThis version of ReadGnssLogger supports Android N\n')
    fprintf('WARNING: did not find "Platform" type in log file, expected "Platform: N"\n')
    fprintf('Please Update GnssLogger\n')
end

v1 = [1;4;0;0];
sCompare = CompareVersions(v,v1);
%Note, we need to check both the logger version (e.g. v1.0.0.0) and the
%Platform version "e.g. Platform: N" for any logic based on version
if strcmp(sCompare,'before')
    fprintf('\nThis version of ReadGnssLogger supports v1.4.0.0 onwards\n')
    error('Found "%s" in log file',line)
end

%% write csv file with header and numbers
gnssCsvfileID = fopen(gnssCsvFilePath,'w');
accCsvfileID = fopen(accCsvFilePath,'w');
gyrCsvfileID = fopen(gyrCsvFilePath,'w');
magCsvfileID = fopen(magCsvFilePath,'w');
nGnss = 0; nAcc = 0; nGyr = 0; nMag = 0;
while ischar(line)
    if contains(line,'Raw,')
        %Now 'line' contains the raw measurements header or data
        line = strrep(line,'Raw,','');
        line = strrep(line,'#',''); line = strrep(line,' ','');
        fprintf(gnssCsvfileID,'%s\n',line);
        nGnss = nGnss + 1;
    elseif contains(line,'UncalAccel,') || contains(line,'Accel,')
        %Now 'line' contains the accel measurements header or data
        if contains(line,'UncalAccel,'), line = strrep(line,'UncalAccel,',''); end
        if contains(line,'Accel,'), line = strrep(line,'Accel,',''); end
        line = strrep(line,'#',''); line = strrep(line,' ','');
        fprintf(accCsvfileID,'%s\n',line);
        nAcc = nAcc + 1;
    elseif contains(line,'UncalGyro,') || contains(line,'Gyro,')
        %Now 'line' contains the gyro measurements header or data
        if contains(line,'UncalGyro,'), line = strrep(line,'UncalGyro,',''); end
        if contains(line,'Gyro,'), line = strrep(line,'Gyro,',''); end
        line = strrep(line,'#',''); line = strrep(line,' ','');
        fprintf(gyrCsvfileID,'%s\n',line);
        nGyr = nGyr + 1;
    elseif contains(line,'UncalMag,') || contains(line,'Mag,')
        %Now 'line' contains the magneto measurements header or data
        if contains(line,'UncalMag,'), line = strrep(line,'UncalMag,',''); end
        if contains(line,'Mag,'), line = strrep(line,'Mag,',''); end
        line = strrep(line,'#',''); line = strrep(line,' ','');
        fprintf(magCsvfileID,'%s\n',line);
        nMag = nMag + 1;
    end
    line = fgetl(txtfileID); % Get next line
end
fclose(txtfileID);
fclose(gnssCsvfileID);
fclose(accCsvfileID);
fclose(gyrCsvfileID);
fclose(magCsvfileID);
if isempty(line) %line should be -1 at eof
    error('\nError occurred while reading file %s\n',fullLogFilePath)
end
% Delete files when only 1 line -> header
if nGnss <= 1
    warning('\t%s does not contain GNSS measurements', fullLogFilePath);
    delete(gnssCsvFilePath)
end
if nAcc <= 1
    warning('\t%s does not contain Accelerometer measurements', fullLogFilePath);
    delete(accCsvFilePath)
end
if nGyr <= 1
    warning('\t%s does not contain Gyroscope measurements', fullLogFilePath);
    delete(gyrCsvFilePath)
end
if nMag <= 1
    warning('\t%s does not contain Magnetometer measurements', fullLogFilePath);
    delete(magCsvFilePath)
end
end %end of function MakeCsv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataTable] = readRawCsv(rawCsvFile)
%% read data from csv file into a numerical matrix 'S' and cell array 'header'
if isfile(rawCsvFile)
    % Columns that are characters
    headerChars = { 'CodeType' };
    % Columns that are longs (int64)
    headerLongs = { 'TimeNanos'
        'FullBiasNanos'
        'ReceivedSvTimeNanos'
        'ReceivedSvTimeUncertaintyNanos'
        'CarrierCycles'
        'elapsedRealtimeNanos'};
    % Import options
    opts = delimitedTextImportOptions('VariableNamesLine', 1, 'DataLines', 2);
    % Preview to obtain column names
    P = preview(rawCsvFile, opts);
    opts.VariableNames = P.Properties.VariableNames;
    clear P
    % Set default type to double
    opts = setvaropts(opts, 'Type', 'double', 'NumberSystem', 'decimal');
    % Change type of char variables if present
    for iVar = 1:length(headerChars)
        if any(contains(opts.VariableNames, headerChars{iVar}))
            opts = setvartype(opts,headerChars{iVar},{'char'});
        end
    end
    % Change type of long variables if present
    for iVar = 1:length(headerLongs)
        if any(contains(opts.VariableNames, headerLongs{iVar}))
            opts = setvartype(opts,headerLongs{iVar},{'int64'});
        end
    end
    
    % Use table to read csv w/o having problems due to chars
    dataTable = readtable(rawCsvFile, opts);
    % Convert char columns into doubles
    auxVec = nan(size(dataTable, 1), length(headerChars));
    for jCol = 1:length(headerChars)
        if any(contains(dataTable.Properties.VariableNames, headerChars{jCol}))
            for iRow = 1:size(dataTable, 1)
                if ~isempty(dataTable.(headerChars{jCol}){iRow})
                    auxVec(iRow, jCol) = double(dataTable.(headerChars{jCol}){iRow});
                end
            end
            dataTable.(headerChars{jCol}) = auxVec(:, jCol);
        end
    end
    clear auxVec
else
    dataTable = [];
end
end% of function ReadRawCsv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gnssRaw,missing] = PackGnssRaw(dataTable)
%% pack data into gnssRaw, and report missing fields
%report clock fields present/missing, based on:
gnssClockFields = {...
    'TimeNanos'
    'TimeUncertaintyNanos'
    'LeapSecond'
    'FullBiasNanos'
    'BiasUncertaintyNanos'
    'DriftNanosPerSecond'
    'DriftUncertaintyNanosPerSecond'
    'HardwareClockDiscontinuityCount'
    'BiasNanos'
    };
%report measurements fields present/missing, based on:
gnssMeasurementFields = {...
    'Cn0DbHz'
    'ConstellationType'
    'MultipathIndicator'
    'PseudorangeRateMetersPerSecond'
    'PseudorangeRateUncertaintyMetersPerSecond'
    'ReceivedSvTimeNanos'
    'ReceivedSvTimeUncertaintyNanos'
    'State'
    'Svid'
    'AccumulatedDeltaRangeMeters'
    'AccumulatedDeltaRangeUncertaintyMeters'
    'BasebandCn0DbHz'
    'FullInterSignalBiasNanos'
    'FullInterSignalBiasUncertaintyNanos'
    'SatelliteInterSignalBiasNanos'
    'SatelliteInterSignalBiasUncertaintyNanos'
    'CodeType'
    'ChipsetElapsedRealtimeNanos'
    };
if ~isempty(dataTable)
    header = dataTable.Properties.VariableNames;
    
    gnssRaw = [];
    
    missing.ClockFields = {};
    missing.MeasurementFields = {};
    
    %pack data into vector variables, if the fields are not NaNs
    for j = 1:length(header)
        if any(isfinite(dataTable{:, j})) %not all NaNs
            gnssRaw.(header{j}) = dataTable{:, j};
        elseif any(strcmp(header{j},gnssClockFields))
            missing.ClockFields{end+1} = header{j};
        elseif any(strcmp(header{j},gnssMeasurementFields))
            missing.MeasurementFields{end+1} = header{j};
        end
    end
else
    fprintf(2,'\tThere are no GNSS measurements available in this file\n')
    gnssRaw = [];
    missing.ClockFields = {'ALL'};
    missing.MeasurementFields = {'ALL'};
end
end %end of function PackGnssRaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meas, missing] = PackMeas(dataTable, type)
%% pack data into meas, and report missing fields
switch type
    case 'acc'
        measFields = {  'utcTimeMillis'
            'elapsedRealtimeNanos'
            'UncalAccelXMps2'
            'UncalAccelYMps2'
            'UncalAccelZMps2'
            'BiasXMps2'
            'BiasYMps2'
            'BiasZMps2'};
    case 'gyr'
        measFields = {  'utcTimeMillis'
            'elapsedRealtimeNanos'
            'UncalGyroXRadPerSec'
            'UncalGyroYRadPerSec'
            'UncalGyroZRadPerSec'
            'DriftXRadPerSec'
            'DriftYRadPerSec'
            'DriftZRadPerSec'};
    case 'mag'
        measFields = {  'utcTimeMillis'
            'elapsedRealtimeNanos'
            'UncalMagXMicroT'
            'UncalMagYMicroT'
            'UncalMagZMicroT'
            'BiasXMicroT'
            'BiasYMicroT'
            'BiasZMicroT'};
end
if ~isempty(dataTable)
    header = dataTable.Properties.VariableNames;
    
    meas = [];
    missing = {};
    
    %pack data into vector variables, if the fields are not NaNs
    for j = 1:length(header)
        if any(isfinite(dataTable{:, j})) %not all NaNs
            meas.(header{j}) = dataTable{:, j};
        elseif any(strcmp(header{j},measFields))
            missing{end+1} = header{j};
        end
    end
else
    fprintf(2, '\tThere are no "%s" measurements available in this file\n', type)
    meas = [];
    missing = {'ALL'};
end
end%end of function PackMeas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gnssRaw,gnssAnalysis,bOk] = CheckGnssClock(gnssRaw,gnssAnalysis)
%% check clock values in gnssRaw
bOk = true;
sFail = ''; %initialize string to record failure messafes
N = length(gnssRaw.ReceivedSvTimeNanos);
%Insist on the presence of TimeNanos (time from hw clock)
if ~isfield(gnssRaw,'TimeNanos')
    s = ' TimeNanos  missing from GnssLogger File.';
    fprintf('WARNING: %s\n',s);
    sFail = [sFail,s];
    bOk = false;
end
if ~isfield(gnssRaw,'FullBiasNanos')
    s = 'FullBiasNanos missing from GnssLogger file.';
    fprintf('WARNING: %s, we need it to get the week number\n',s);
    sFail = [sFail,s];
    bOk = false;
end
if ~isfield(gnssRaw,'BiasNanos')
    gnssRaw.BiasNanos = zeros(N,1);
end
if ~isfield(gnssRaw,'HardwareClockDiscontinuityCount')
    %possibly non fatal error, we assume there is no hardware clock discontinuity
    %so we set to zero and move on, but we print a warning
    gnssRaw.HardwareClockDiscontinuityCount = zeros(N,1);
    fprintf('WARNING: Added HardwareClockDiscontinuityCount=0 because it is missing from GNSS Logger file\n');
end

%check FullBiasNanos, it should be negative values
bChangeSign = any(gnssRaw.FullBiasNanos<0) & any(gnssRaw.FullBiasNanos>0);
assert(~bChangeSign,...
    'FullBiasNanos changes sign within log file, this should never happen');
%Now we know FullBiasNanos doesnt change sign,auto-detect sign of FullBiasNanos,
%if it is positive, give warning and change
if any(gnssRaw.FullBiasNanos>0)
    gnssRaw.FullBiasNanos = -1*gnssRaw.FullBiasNanos;
    fprintf('WARNING: FullBiasNanos wrong sign. Should be negative. Auto changing inside ReadGpsLogger\n');
    gnssAnalysis.GnssClockErrors = [gnssAnalysis.GnssClockErrors,...
        sprintf(' FullBiasNanos wrong sign.')];
end

%compute full cycle time of measurement, in milliseonds
gnssRaw.allRxMillis = int64((gnssRaw.TimeNanos - gnssRaw.FullBiasNanos)*1e-6);
%allRxMillis is now accurate to one millisecond (because it's an integer)

if ~bOk
    gnssAnalysis.ApiPassFail = ['FAIL ',sFail];
end
end %end of function CheckGnssClock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GnssAnalysis, SensorAnalysis] = ReportMissingFields(GnssAnalysis,missing, accMissing, gyrMissing, magMissing)
%% report missing clock and measurement fields in gnssAnalysis

%report missing clock fields
if ~isempty(missing.ClockFields)
    GnssAnalysis.GnssClockErrors = sprintf(...
        '%s Missing Fields:',GnssAnalysis.GnssClockErrors);
    for i=1:length(missing.ClockFields)
        GnssAnalysis.GnssClockErrors = sprintf(...
            '%s %s,',GnssAnalysis.GnssClockErrors,missing.ClockFields{i});
    end
    GnssAnalysis.GnssClockErrors(end) = '.';%replace final comma with period
end

%report missing measurement fields
if ~isempty(missing.MeasurementFields)
    GnssAnalysis.GnssMeasurementErrors = sprintf(...
        '%s Missing Fields:',GnssAnalysis.GnssMeasurementErrors);
    for i=1:length(missing.MeasurementFields)
        GnssAnalysis.GnssMeasurementErrors = sprintf(...
            '%s %s, ',GnssAnalysis.GnssMeasurementErrors,...
            missing.MeasurementFields{i});
    end
    GnssAnalysis.GnssMeasurementErrors(end-1:end) = '. ';%replace last comma with period
end

%assign pass/fail
if ~any(strfind(GnssAnalysis.ApiPassFail,'FAIL')) %didn't already fail
    if isempty(missing.ClockFields) && isempty(missing.MeasurementFields)
        GnssAnalysis.ApiPassFail = 'PASS';
    else
        GnssAnalysis.ApiPassFail = 'FAIL BECAUSE OF MISSING FIELDS';
    end
end

%report missing sensors' fields
SensorAnalysis.accErrors = '';
SensorAnalysis.gyrErrors = '';
SensorAnalysis.magErrors = '';
if ~isempty(accMissing)
    SensorAnalysis.accErrors = 'Accelerometer Errors. Missing Fields: ';
    for i=1:length(accMissing)
        SensorAnalysis.accErrors = [SensorAnalysis.accErrors accMissing{i} ', '];
    end
    SensorAnalysis.accErrors(end-1:end) = '. ';%replace last comma with period
end
if ~isempty(gyrMissing)
    SensorAnalysis.gyrErrors = 'Gyroscope Errors. Missing Fields: ';
    for i=1:length(gyrMissing)
        SensorAnalysis.gyrErrors = [SensorAnalysis.gyrErrors gyrMissing{i} ', '];
    end
    SensorAnalysis.gyrErrors(end-1:end) = '. ';%replace last comma with period
end
if ~isempty(magMissing)
    SensorAnalysis.magErrors = 'Magnetometer Errors. Missing Fields: ';
    for i=1:length(magMissing)
        SensorAnalysis.magErrors = [SensorAnalysis.magErrors magMissing{i} ', '];
    end
    SensorAnalysis.magErrors(end-1:end) = '. ';%replace last comma with period
end
end %end of function ReportMissingFields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2016 Google Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


