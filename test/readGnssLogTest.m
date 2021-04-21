clear; close all; clc;

%% Configuration
% Measurement campaign name
campaignName = '2020-08-06-US-MTV-2'; %'2020-07-08-US-MTV-1';
% Phone name
phoneName = 'Pixel4XL'; 
% (num of Raw) MINUS ONE in file: grep -o -i Raw Pixel4XL_GnssLog.txt | wc -l
numRaw = 49407; %'2020-07-08-US-MTV-1', Pixel4Xl: 67818; 

%% Dataset
% Dataset path
datasetsPath = './data/training/datasets/';
testFileName = [phoneName '_GnssLog.txt'];
derivedFileName = [phoneName '_GnssLog.derived'];
dirName = [datasetsPath campaignName '/'];

% %% Read data with FRANK's code
% [gnssRawRef, gnssAnalysisRef] = ReadGnssLogger(dirName,testFileName);

%% Read data with NEW code
[gnssRawNew, nssAnalysisNew, accRawNew, gyrRawNew, magRawNew, sensorAnalysis] = readGnssLog(dirName,testFileName);

% %% Test 1: Check if new gnssRaw structure is equal to the one obtained with Frank's code
% if isequaln(gnssRawRef, gnssRawNew), fprintf(' ======== Test 1 - OK ========');
% else, error('   Test 1 - Failed. New gnssRaw differs from Ref gnssRaw');
% end

%% Test 1: Check correct num of measurements
if length(gnssRawNew.utcTimeMillis) == numRaw, disp(' ======== Test 1 - OK ========');
else, error('   Test 1 - Failed. Size is not correct');
end

%% Test 2: Check accelerometer values
accFields = {   'utcTimeMillis'
                'elapsedRealtimeNanos'
                'UncalAccelXMps2'
                'UncalAccelYMps2'
                'UncalAccelZMps2'
                'BiasXMps2'
                'BiasYMps2'
                'BiasZMps2'};
accValsFirst = [];
accValsLast = [];
%'2020-07-08-US-MTV-1', Pixel4Xl:
% accValsFirst = [1594247262456,17373972330561,0.49355114,9.520887,0.86991197,0.0,0.0,0.0];
% accValsLast = [1594249128009,19239526563068,-0.88991547,9.36713,0.101608306,0.0,0.0,0.0];

if isempty(accValsFirst) && isempty(accRawNew)
    disp(' ======== Test 2 - OK ========');
elseif ~isempty(accValsFirst) && ~isempty(accRawNew)
    % Test 2.1: Check existence of fields
    isOk = 1;
    for i = 1:length(accFields), isOk = isOk*isfield(accRawNew, accFields{i}); end

    if isOk, disp(' ======== Test 2.1 - OK ========');
    else, error('   Test 2.1 - Failed. accRawNew does not contain all the fields');
    end
    clear isOk;

    % Test 2.2: Check first and last set of measurements
    isFirstOk = 1; isLastOk = 1;
    for i = 1:length(accFields), isFirstOk = isFirstOk*(accRawNew.(accFields{i})(1) == accValsFirst(i)); end
    for i = 1:length(accFields), isLastOk = isLastOk*(accRawNew.(accFields{i})(end) == accValsLast(i)); end

    if isFirstOk && isLastOk,   disp(' ======== Test 2.2 - OK ========');
    elseif isFirstOk,           error('   Test 2.2 - Failed. Last set of measurements is wrong');
    elseif isLastOk,            error('   Test 2.2 - Failed. First set of measurements is wrong');
    else,                       error('   Test 2.2 - Failed. First and last set of measurements are wrong');
    end
    clear isFirstOk; clear isLastOk;
else
    error('   Test 2 - Failed. Dimensions do not match between reference and obtained measurements');
end

%% Test 3: Check gyro values
gyrFields = {   'utcTimeMillis'
                'elapsedRealtimeNanos'
                'UncalGyroXRadPerSec'
                'UncalGyroYRadPerSec'
                'UncalGyroZRadPerSec'
                'DriftXRadPerSec'
                'DriftYRadPerSec'
                'DriftZRadPerSec'};
            
gyrValsFirst = [];
gyrValsLast = [];

%'2020-07-08-US-MTV-1', Pixel4Xl:
% gyrValsFirst = [1594247262455,17373971210561,-0.0815042,-0.036781326,0.045842,-3.9329554E-4,-2.0115673E-4,0.0056578554];
% gyrValsLast = [1594249127999,19239515755412,-0.059221078,0.11180529,0.04539248,-9.5147477E-4,3.4031557E-4,0.0066613476];

if isempty(gyrValsFirst) && isempty(gyrRawNew)
    disp(' ======== Test 3 - OK ========');
elseif ~isempty(gyrValsFirst) && ~isempty(gyrRawNew)
    % Test 3.1: Check existence of fields
    isOk = 1;
    for i = 1:length(gyrFields), isOk = isOk*isfield(gyrRawNew, gyrFields{i}); end

    if isOk, disp(' ======== Test 3.1 - OK ========');
    else, error('   Test 3.1 - Failed. gyrRawNew does not contain all the fields');
    end
    clear isOk;

    % Test 3.2: Check first and last set of measurements
    isFirstOk = 1; isLastOk = 1;
    for i = 1:length(gyrFields), isFirstOk = isFirstOk*(gyrRawNew.(gyrFields{i})(1) == gyrValsFirst(i)); end
    for i = 1:length(gyrFields), isLastOk = isLastOk*(gyrRawNew.(gyrFields{i})(end) == gyrValsLast(i)); end

    if isFirstOk && isLastOk,   disp(' ======== Test 3.2 - OK ========');
    elseif isFirstOk,           error('   Test 3.2 - Failed. Last set of measurements is wrong');
    elseif isLastOk,            error('   Test 3.2 - Failed. First set of measurements is wrong');
    else,                       error('   Test 3.2 - Failed. First and last set of measurements are wrong');
    end
    clear isFirstOk; clear isLastOk;
else
    error('   Test 3 - Failed. Dimensions do not match between reference and obtained measurements');
end

%% Test 4: Check magneto values
magFields = {   'utcTimeMillis'
                'elapsedRealtimeNanos'
                'UncalMagXMicroT'
                'UncalMagYMicroT'
                'UncalMagZMicroT'
                'BiasXMicroT'
                'BiasYMicroT'
                'BiasZMicroT'};

magValsFirst = [];
magValsLast = [];

%'2020-07-08-US-MTV-1', Pixel4Xl:
% magValsFirst = [1594247262455,17373971618738,5.4223285,-84.088234,-101.32536,-5.7357035,-51.29586,-91.78012];
% magValsLast = [1594249128012,19239528908589,23.372301,-87.13828,-71.227104,-5.7357035,-51.29586,-91.78012];

if isempty(magValsFirst) && isempty(magRawNew)
    disp(' ======== Test 4 - OK ========');
elseif ~isempty(magValsFirst) && ~isempty(magRawNew)
    % Test 4.1: Check existence of fields
    isOk = 1;
    for i = 1:length(magFields), isOk = isOk*isfield(magRawNew, magFields{i}); end

    if isOk, disp(' ======== Test 4.1 - OK ========');
    else, error('   Test 4.1 - Failed. magRawNew does not contain all the fields');
    end
    clear isOk;

    % Test 4.2: Check first and last set of measurements
    isFirstOk = 1; isLastOk = 1;
    for i = 1:length(magFields), isFirstOk = isFirstOk*(magRawNew.(magFields{i})(1) == magValsFirst(i)); end
    for i = 1:length(magFields), isLastOk = isLastOk*(magRawNew.(magFields{i})(end) == magValsLast(i)); end

    if isFirstOk && isLastOk,   disp(' ======== Test 4.2 - OK ========');
    elseif isFirstOk,           error('   Test 4.2 - Failed. Last set of measurements is wrong');
    elseif isLastOk,            error('   Test 4.2 - Failed. First set of measurements is wrong');
    else,                       error('   Test 4.2 - Failed. First and last set of measurements are wrong');
    end
    clear isFirstOk; clear isLastOk;
else
    error('   Test 4 - Failed. Dimensions do not match between reference and obtained measurements');
end
