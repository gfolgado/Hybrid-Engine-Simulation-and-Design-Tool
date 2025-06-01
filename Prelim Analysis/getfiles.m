function [FileList, Pc] = getfiles()
%GETFILES Returns array with list of .out files and Pc values

cd('CEA');

files = dir;
exclusions = {'.', '..', 'cea.inc', 'cea2.csv', 'cea2.f', 'cea2.inp', ...
    'cea2.out', 'cea2.plt', 'cea2m.f', 'cea2m.inc', 'FCEA2.exe', ...
    'FCEA2m.exe', 'Guide Analysis.pdf', 'test.csv', 'thermo.inp', ...
    'thermo.lib', 'thermo.out', 'trans.inp', 'trans.lib', 'trans.out', ...
    'Users Manual.pdf', 'runCEA.bat'};

FileList = {};
Pc = [];

for i = 1:numel(files)
    if ~files(i).isdir && ~any(strcmp(files(i).name, exclusions))
        if ~endsWith(files(i).name, '.inp')
            current = files(i).name;
            current = erase(current, '.out');
            FileList = [FileList, current];
            press_str = extractBetween(current, 'press', '.out');
            press_val = str2double(press_str);
            if ~isnan(press_val)
                Pc = [Pc, press_val];
            end
        end
    end
end

Pc = sort(Pc);
FileList = sort(FileList);

cd('..');
end

