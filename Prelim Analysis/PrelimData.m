function prelim_data = PrelimData(fileNames)
%PrelimData Collects data from evey .out file
%   Collects performance data for every generated .out file

prelim_data =[];

outs = strcat(fileNames, '.out');

for i = 1:numel(outs)
    
   performance_matrix = CEAPrelimMatrix(outs{i});
   prelim_data(:,:,i) = performance_matrix;
    
end
end

