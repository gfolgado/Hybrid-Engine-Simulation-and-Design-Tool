function error = AutoRunCEA(fileNames)
%UNTITLED Creates a .bat file to run CEA for every file in fileNames
%   Returns 1 in case of error detected. Creates a CEA .out file for every
%   .inp file in the fileNames cell array

error = 0;
% Folder with CEA and .inp files
folderName = 'CEA';

% Batch file to be created to automate CEA runs
batch = 'runCEA.bat';

% Open CEA directory
cd(folderName);

try
    % Creates batch file, returns error=1 if file cannot be created
    fid = fopen(batch, 'w');
    if fid == -1
        error = 1;
        return
    end

    
    % Goes through all file names in array
    for i = 1: numel(fileNames)
        % Creates line of code in batch file to run CEA .exe and use file
        % name to indicate which calculation to perform
        fprintf(fid,'echo %s | FCEA2m.exe\n', fileNames{i});
    end
    
    % Closes batch file
    fclose(fid);
    % Runs batch file
    system('runCEA.bat > nul 2>&1');
    
    
catch
    % Makes error=1 if some issue occurs in code    
    error = 1;
end

% Returns to tool directory
cd('..');
clc
end

