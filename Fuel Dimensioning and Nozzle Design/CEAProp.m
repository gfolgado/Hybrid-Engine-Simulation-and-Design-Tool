function data = CEAProp(Pc, OF, eps, fuels, oxs)
%CEAPROP Runs CEA for input values, organizes results into vector for
%HybridPropulsion function use

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

% Runs CEAMatrixv2, to gather data from the
% CEA output files. Then this data is organized in a vector pertinent to
% the HybridPropulsion.m function

data_full = CEAMatrixv2(FileList(1));
data = [data_full(6), data_full(8), data_full(10), data_full(12), data_full(9)];
end

