function Derived = readDerivedData(dirName, fileName)
% READDERIVEDDATA reads the .derived data files that can be obtained from
% GNSS Analysis tool and returns the data in a structure.
%
%   Derived = READDERIVEDDATA(dirName, fileName)
% 
% Input:  
%   dirName = name of the directory.
%   fileName = name of the file to read.
% Output:  
%   Derived = structure containing the derived data. 
%
% The Different types of data are saved as tables in different attributes. 
% The names of the attributes are set the same as in the data types in the 
% input file.

% Author: Arnau Ochoa Banuelos (CS Group), March 2021

%% Initializations
% Lines containing the data description
dataDescriptionLines = [6 11];

%% Read header
derivedHeaderOpts = delimitedTextImportOptions();
derivedHeaderOpts.Whitespace = '#';                 % Convert # into whitespace
derivedHeaderOpts.DataLines = dataDescriptionLines; 
derivedHeaderTable = readtable([dirName fileName], derivedHeaderOpts);
numCols = size(derivedHeaderTable, 2);

%% Read data
derivedDataOpts = delimitedTextImportOptions('NumVariables', numCols);
derivedDataOpts.CommentStyle = '#';                                 % Ignore comments
derivedDataOpts = setvartype(derivedDataOpts, 2:numCols, 'double'); % Set all columns but the first as doubles
derivedDataTable = readtable([dirName fileName], derivedDataOpts);


%% Refactor header
derivedHeader = cell(size(derivedHeaderTable, 1), 1);
for iRow = 1:size(derivedHeaderTable, 1)
    % Remove whitespaces at the beginning of the header
    derivedHeaderTable{iRow,1}{1}(1) = [];
    % Transform to cell array
    derivedHeader{iRow} = derivedHeaderTable{iRow, :};
    % Remove empty cells
    derivedHeader{iRow, :} = derivedHeader{iRow, :}(~cellfun('isempty',derivedHeader{iRow, :}));
end

%% Refactor data
Derived = struct;
for iRow = 1:length(derivedHeader)
    % Put table into new struct field
    Derived.(derivedHeader{iRow}{1}) = derivedDataTable(strcmp(derivedDataTable.Var1, derivedHeader{iRow}{1}), :);
    % Remove empty columns
    Derived.(derivedHeader{iRow}{1}) = Derived.(derivedHeader{iRow}{1})(:, 2:length(derivedHeader{iRow}));
    % Assign correct names from header description
    Derived.(derivedHeader{iRow}{1}).Properties.VariableNames = derivedHeader{iRow}(2:end);
end

end