function [tables, FileList] = CEA_TableGenerator(Pc, OF, eps, fuels, oxs)
%CEA_TABLEGENERATOR Generates tables for thermo data for CC model
%   A 3D matrix is created, where each 2D matrix contains both thermal and
%   performance data from CEA for a specific Pc value, in a range of OF
%   ratio values. Check the header comments for CEAMatrix.m function for
%   the index of the 2D matrices.

% Creates .inp files for automatic running of CEA
[error, FileList] = CEAMakeInp (Pc,OF,eps,fuels,oxs);

if error == 1
    fprintf('An error was detected during the .inp file creation process');
    return;
end

% Creates and runs batch file to run all .inp files created
error = AutoRunCEA(FileList);

if error == 1
    fprintf('An error was detected during the batch file creation/execution')
    return;
end

% For every .inp file (meaning, for every chamber pressure value in Pc),
% the function CEAMatrix is executed, and the resulting 2D matrix is added
% to a 3D matrix, which is the output of this function

tables = [];

    for i = 1:numel(FileList)

        table = CEAMatrix(FileList(i));
        tables(:,:,i) = table;

    end
end

