function [FileList, prelim_data] = CEAPrelimMaster(Pc,OF,PAtm,fuels, oxs)
%CEAPrelimMaster Takes user inputs into CEA to perform prelim analysis
%   Calls function CEAMakeInp to create input files for CEA analysis, then
%   calls AutoRunCEA to create a .bat file to run CEA for all those created
%   .inp files, and finally PrelimData function is called to get the 3D
%   matrix of CEA results, for easy plotting of performance parameters.

[error, FileList] = CEAMakeInp(Pc, OF, PAtm, fuels, oxs);
if error == 1
   fprintf('An error ocurred while creating .inp files\n');
   return
end

error = AutoRunCEA(FileList);

if error == 1
   fprintf('An error ocurred while running CEA\n');
   return
end

prelim_data = PrelimData(FileList);
end

