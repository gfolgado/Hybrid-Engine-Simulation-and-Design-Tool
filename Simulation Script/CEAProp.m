function [dataProp, dataNozzle] = CEAProp(Pc, OF, eps, fuels, oxs)
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
% the HybridPropulsion.m function (dataProp, which holds in order values 
% for T_c, rho_c, cp_c, k_c and MW_c), and another vector named dataNozzle
% which holds data pertinent to the nozzle flow model not present in dataProp
% (Me_sub, Me_sup, pe_sub, pe_sup, cstar_sup)

data_full = CEAMatrixv2(FileList(1));
dataProp = [data_full(6), data_full(8), data_full(10), data_full(12), data_full(9)];
dataNozzle = [data_full(13), data_full(14), data_full(4), data_full(5), data_full(16)];
end

