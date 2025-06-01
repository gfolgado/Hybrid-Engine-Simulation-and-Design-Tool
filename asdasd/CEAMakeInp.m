function [error, FileList] = CEAMakeInp(Pc, OF, Patm, fuels, oxs)
%   CEAMakeInp creates .inp files for running CEA Software analysis
%   The user gives as inputs: array of chamber pressure values (Pc); an
%   array of OF ratios he wishes to analyse; the value of Patm for which
%   the nozzle will be perfectly expanded; fuels array of cells where the
%   first row is the name of the reagents, the second the weight fraction 
%   of each species, the temperature of the fuels, and if the reagents are 
%   not present in the CEA reagents library, the forth is the specific 
%   enthalpy, and the fifth the exploded chemical formula of the reagents
%   (it should be noted that, for chemical species in the CEA library, the 
%   value of enthalpy is set to zero, and the exploded formula is an empty
%   string); finally an array of cells for the oxidizers, arranged in
%   the same fashion as the fuels array. The output of this function is the
%   creating of .inp files for CEA (one for each chamber pressure value in 
%   the array), an error flag, which will be 0 if no error is detected, and
%   1 if there is, and a list of all the file names created when running
%   the function

error = 0;
folder = 'CEA';
FileList = {};

% If there is no folder named CEA, one is created
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Check if there are files in the folder, deletes them
files = dir(fullfile(folder, '*'));
% Exclude all CEA files so they are not deleted when clearing old .inp files
excl = {'.', '..', 'cea.inc', 'cea2.csv', 'cea2.f', 'cea2.inp', ...
    'cea2.out', 'cea2.plt', 'cea2m.f', 'cea2m.inc', 'FCEA2.exe', ...
    'FCEA2m.exe', 'Guide Analysis.pdf', 'test.csv', 'thermo.inp', ...
    'thermo.lib', 'thermo.out', 'trans.inp', 'trans.lib', 'trans.out', ...
    'Users Manual.pdf'};

if ~isempty(files)
    fileNames = {files.name};
    for j = 1:numel(fileNames)
        % Initialize flag to indicate whether the current file is excluded
        excluded = false;
        % Check if the current file name matches any of the excluded files
        for k = 1:numel(excl)
            if strcmp(fileNames{j}, excl{k})
                % Set flag to true if match found
                excluded = true;
                break; % No need to continue checking if match is found
            end
        end
        % If the current file is not in the list of excluded files, delete it
        if ~excluded
            delete(fullfile(folder, fileNames{j}));
        end
    end
end

for i = 1:length(Pc)

    % Generate the file name
    file = fullfile(folder, ['press', num2str(Pc(i)), '.inp']);
    name = ['press', num2str(Pc(i))];
    FileList(i) = cellstr(name); 

    % Create an empty .inp file
    fid = fopen(file, 'w');
    
    % Write in .inp file
    fprintf(fid, 'problem\no/f=');
    
    % Length of the line after writing o/f=
    linelength = 4;
    
    % Go through every OF value in array to write in file
    for j = 1:length(OF)
        
        % Calculates line length IF new value is printed
        linelength = linelength + numel(sprintf('%.4f', OF(j))) + 1; 
        
        % If line length will go over 40, break line, updates linelength to
        % be the length of the OF value that will be written
        if linelength > 40
            fprintf(fid, '\n');
            linelength = numel(sprintf('%.4f', OF(j)));
        end 
        
        % Writes new OF value
        fprintf(fid, '%.4f,', OF(j));        
    end
    
    % Write problem description, only rocket problems with frozen flow are
    % considered for this preliminary application
    fprintf(fid,'\n    rocket  eql\n');
    fprintf(fid,'  p,bar=%.4f,\n', Pc(i));
    
    % Calculate inverse pressure ratio to use as input in CEA (enables
    % preliminary calculation of Ae/At for optimum expansion at given Patm)
    pip = Pc(i) / Patm;
    
    % Writes computed inverse pressure ratio in .inp file
    fprintf(fid,'  pi/p=%.4f,\n', pip);
    
    % Introduce reagents inputs
    fprintf(fid,'react\n');
    
    
    % Checks if fuel weight fractions sum is equal to one, otherwise, 
    % return error as 1
    weights = cell2mat(fuels(:, 2));
    if sum(weights) ~= 1
        error = 1;
    end
    
    % Same thing for oxidizers
    weights = cell2mat(oxs(:, 2));
    if sum(weights) ~= 1
        error = 1;
    end
    
    % Checks if temperature values are not negative, since values are in K
    for j=1:size(fuels, 1)
        if fuels{j,3} < 0
            error = 1;
        end
    end
    
    for j=1:size(oxs, 1)
        if oxs{j,3} < 0
            error = 1;
        end
    end
    
    % Go through fuel information matrix, write info in .inp file
    for j = 1:size(fuels, 1)
        species = fuels{j, 1};
        weight_fraction = fuels{j, 2};
        temp = fuels{j, 3};

        % Check if exploded formula is empty string (reagent in library)
        if isempty(fuels{j, 5})
            % If exploded formula is empty, write only fuel information
            fprintf(fid, '  fuel=%s wt=%.3f t,k=%.3f\n', species, weight_fraction, temp);
        else
            % If exploded formula is not empty, write additional information
            h_specific = fuels{j, 4};
            exploded_formula = fuels{j, 5};
            fprintf(fid, '  fuel=%s wt=%.3f t,k=%.3f\n  h,kj/mol=%.3f %s\n', species, weight_fraction, temp, h_specific, exploded_formula);
        end
    end
    
    % Same thing for oxidizers list
    for j=1:size(oxs, 1)
        species = oxs{j, 1};
        weight_fraction = oxs{j,2};
        temp = oxs{j,3};
        
        % Check if exploded formula is empty string (reagent in library)
        if isempty(oxs{j, 5})
            % If exploded formula is empty, write only fuel information
            fprintf(fid, '  oxid=%s wt=%.3f t,k=%.3f\n', species, weight_fraction, temp);
        else
            % If exploded formula is not empty, write additional information
            h_specific = oxs{j, 4};
            exploded_formula = oxs{j, 5};
            fprintf(fid, '  oxid=%s wt=%.3f t,k=%.3f\n  h,kj/mol=%.3f %s\n', species, weight_fraction, temp, h_specific, exploded_formula);
        end
    end
    
    fprintf(fid,'output siunits short\nend');
    
    fclose(fid);
end
end

