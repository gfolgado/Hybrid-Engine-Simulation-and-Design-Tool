function [Pc, OF, Patm, fuels, oxs, reagents] = UserInputs()
%USERINPUTS Queries user for inputs to use in CEA prelim analysis

%% Atmospheric Model

% To calculate inverse pressure ratios, to obtain preliminary values of
% optimum nozzle expansion ratio, the atmospheric pressure at the altitude
% the user desires to have the nozzle perfectly expanded must be computed

clc;
h_site_ASL = input('Enter launch site altitude ASL [m]: '); 
h_burnout = input('Enter estimated altitude AGL of burnout [m]: ');
Pc = [];
OF = [];


h_opt = h_site_ASL + h_burnout / 2;
p_a_site = (101325 * (1 - 2.25577 * 10^-5 * h_site_ASL)^5.25588) * 10^-5;
Patm = (101325 * (1 - 2.25577 * 10^-5 * h_opt)^5.25588) * 10^-5;


%% User Inputs (Chamber Pressures and OF ratios)

while true
    % Prompt the user to choose the input method
    input_method = input('Do you wish to input values for chamber pressure one at a time (type "1") or indicate a range (type "2")? ', 's');
    clc

    if strcmp(input_method, '1')  % Input values one at a time
        while true
            % Prompt the user for a chamber pressure value
            input_str = input('Enter chamber pressure values (one at a time, must be numerical; write a negative number to stop): ', 's');
            clc

            % Convert input string to number (returns NaN if not numerical)
            current = str2double(input_str);

            % Check if the input is not a valid numerical value (NaN)
            if isnan(current)
                % Display an error message for invalid input
                fprintf('Invalid input. Please enter a numerical value.\n');
            % Check if the input is negative or zero (signaling the end of input)
            elseif current <= 0
                break;  % Exit the loop if input is negative or zero
            else
                % Append the input value to the array of chamber pressures
                Pc = [Pc, current];
            end
        end
        break;  % Exit the input method loop

    elseif strcmp(input_method, '2')  % Input a range
        while true
            % Prompt the user for the range values
            min_value = input('Enter the minimum value for the chamber pressure range (must be positive): ');
            max_value = input('Enter the maximum value for the chamber pressure range (must be larger than minimum): ');
            num_elements = input('Enter the number of elements for the chamber pressure range (must be integer greater or equal to 2): ');
            clc

            % Validate the input values
            if min_value <= 0
                fprintf('Invalid minimum value. Please enter a positive value.\n');
            elseif max_value <= min_value
                fprintf('Invalid maximum value. Please enter a value greater than the minimum.\n');
            elseif num_elements < 2 || floor(num_elements) ~= num_elements
                fprintf('Invalid number of elements. Please enter a integer equal to or greater than 2.\n');
            else
                % Generate the Pc values using linspace
                Pc = linspace(min_value, max_value, num_elements);
                break;  % Exit the input loop
            end
        end
        break;  % Exit the input method loop

    else
        % Display an error message for invalid input method and continue the loop
        fprintf('Invalid input method selected. Please enter either "1" or "2".\n');
    end
end

clc
Pc = sort(Pc);

while true
    % Prompt the user to choose the input method
    input_method = input('Do you wish to input values for OF ratio one at a time (type "1") or indicate a range (type "2")? ', 's');
    clc

    if strcmp(input_method, '1')  % Input values one at a time
        while true
            % Prompt the user for an OF value
            input_str = input('Enter an OF value (must be numerical; write a negative number to stop): ', 's');
            clc

            % Convert input string to number (returns NaN if not numerical)
            current = str2double(input_str);

            % Check if the input is not a valid numerical value (NaN)
            if isnan(current)
                % Display an error message for invalid input
                fprintf('Invalid input. Please enter a numerical value.\n');
            % Check if the input is negative (signaling the end of input)
            elseif current < 0
                break;  % Exit the loop if input is negative
            else
                % Append the input value to the array of OF ratios
                OF = [OF, current];
            end
        end
        break;  % Exit the input method loop

    elseif strcmp(input_method, '2')  % Input a range
        while true
            % Prompt the user for the range values
            min_value = input('Enter the minimum value for the OF range (must be non-negative): ');
            max_value = input('Enter the maximum value for the OF range (must be greater than minimum value): ');
            clc

            % Validate the input values
            if min_value < 0
                fprintf('Invalid minimum value. Please enter a non-negative value.\n');
            elseif max_value <= min_value
                fprintf('Invalid maximum value. Please enter a value greater than the minimum.\n');
            else
                % Generate the OF values using linspace
                num_elements = 50; % Fixed number of elements for OF array
                OF = linspace(min_value, max_value, num_elements);
                break;  % Exit the input loop
            end
        end
        break;  % Exit the input method loop

    else
        % Display an error message for invalid input method and continue the loop
        fprintf('Invalid input method selected. Please enter either "1" or "2".\n');
    end
end

OF = sort(OF);

%% Get CEA Reagents List

% Initialize reagent array
reagents = {};

% Open CEA directory
cd('CEA');

try
    
    % Open thermo.out file which contains all reagents names
    fid = fopen('thermo.out', 'r');
    
    if fid == -1
        error('Error: Could not open the thermo.out file.');
    end
    
    
    % Skips all lines until 'thermo' string is found
    tline = fgetl(fid);
    while ischar(tline)
        if contains(tline, 'thermo')
            break;  % Exit the loop if 'thermo' is found
        end
       tline = fgetl(fid);  % Read the next line
    end
    
    tline = fgetl(fid); % Moves to first reagent line
    
    % Get reagent names and save them in array
    while ischar(tline)
        
    % Remove leading whitespace
    tline = strtrim(tline);
    
    % Split the line by spaces
    parts = strsplit(tline, ' ');
    
    % Extract the first element (reagent name)
    current = parts{1};
    
    % Append the reagent name to the cell array
    reagents = [reagents, current];
    
    % Read the next line
    tline = fgetl(fid);
    end
    
catch
    % Error message if something unexpected happens
    error('Error: An error occurred while reading the thermo.out file.');
end

% Return to MATLAB code directory
cd('..');

%% Fuels and Oxidizers Input

clc;
flag = 1;

while flag
    fuels = {}; % Reset fuels matrix for each iteration
    
    while true
        % Prompt the user for the fuel name
        fuel_name = input('Enter the fuel name: ', 's'); 

        % Check if the fuel name matches any reagents in the CEA library
        if any(strcmp(fuel_name, reagents))
            % Input weight fraction
            wt = input('Enter fuel weight fraction (note that this value goes from 0 to 1,\nand the sum of all weight fractions if one or more fuels are selected must be equal to 1): ');

            % Redo weight fraction input if invalid input
            while ~(isnumeric(wt) && wt > 0 && wt <= 1)
                clc;
                fprintf('Invalid weight fraction value.\n');
                wt = input('Enter the fuel weight fraction: '); 
            end

            % Input reagent temp
            temp = input('Enter the fuel temperature in K: ');

            % Redo temp input if invalid input
            while ~(isnumeric(temp) && temp > 0)
                clc;
                fprintf('Invalid temperature value.\n');
                temp = input('Enter the fuel temperature in K: '); 
            end
            
            % Add empty cells of enthalpy and exploded formula for array
            % dimension consistency
            
            h_specific = 0;
            exploded_formula = '';
            
            % Add fuel and properties to matrix
            fuels = [fuels; {fuel_name, wt, temp, h_specific, exploded_formula}];
        else
            % Prompt the user if they wish to use the fuel anyway
            prompt = sprintf('Fuel name does not match any reagents in CEA library. Do you wish to use this fuel anyway? (y/n): ');
            user_input = input(prompt, 's');

            if strcmpi(user_input, 'y')
                % Input weight fraction
                wt = input('Enter fuel weight fraction (note that this value goes from 0 to 1,\nand the sum of all weight fractions if one or more fuels are selected must be equal to 1): ');

                % Redo weight fraction input if invalid input
                while ~(isnumeric(wt) && wt > 0 && wt <= 1)
                    clc;
                    fprintf('Invalid weight fraction value.\n');
                    wt = input('Enter the fuel weight fraction: '); 
                end

                % Input reagent temp
                temp = input('Enter the fuel temperature in K: ');

                % Redo temp input if invalid input
                while ~(isnumeric(temp) && temp > 0)
                    clc;
                    fprintf('Invalid temperature value.\n');
                    temp = input('Enter the fuel temperature in K: '); 
                end

                % Input specific enthalpy in Kj/mol
                h_specific = input('Enter the specific enthalpy in Kj/mol: ');

                % Input exploded chemical formula
                exploded_formula = input('Enter the exploded chemical formula:\nExample: if species is CO2, input should be C 1 O 2\n', 's');

                % Add fuel and properties to matrix
                fuels = [fuels; {fuel_name, wt, temp, h_specific, exploded_formula}];
            else
                % If the user chooses not to use the fuel, prompt again for fuel name
                continue;
            end
        end

        % Option to include more fuels
        more = input('Do you wish to add another fuel? (y/n): ', 's');
        if ~strcmpi(more, 'y')
            break;
        end
        clc
    end

    % Calculate the sum of weight fractions
    sum_weight_fractions = sum([fuels{:, 2}]);

    % Check if the sum is not equal to 1
    if sum_weight_fractions ~= 1
        fprintf('Error: Sum of weight fractions is not 1. Please input the fuels again.\n');
    else
        flag = 0; % Redo fuel inputs if sum of weight fractions is not 1
    end

end

flag = 1;
clc;

clc;
flag = 1;

while flag
    oxs = {}; % Reset oxs matrix for each iteration
    
    while true
        % Prompt the user for the ox name
        ox_name = input('Enter the ox name: ', 's'); 

        % Check if the ox name matches any reagents in the CEA library
        if any(strcmp(ox_name, reagents))
            % Input weight fraction
            wt = input('Enter ox weight fraction (note that this value goes from 0 to 1,\nand the sum of all weight fractions if one or more oxs are selected must be equal to 1): ');

            % Redo weight fraction input if invalid input
            while ~(isnumeric(wt) && wt > 0 && wt <= 1)
                clc;
                fprintf('Invalid weight fraction value.\n');
                wt = input('Enter the ox weight fraction: '); 
            end

            % Input reagent temp
            temp = input('Enter the ox temperature in K: ');

            % Redo temp input if invalid input
            while ~(isnumeric(temp) && temp > 0)
                clc;
                fprintf('Invalid temperature value.\n');
                temp = input('Enter the ox temperature in K: '); 
            end
            
            % Add empty cells of enthalpy and exploded formula for array
            % dimension consistency
            
            h_specific = 0;
            exploded_formula = '';
            
            % Add ox and properties to matrix
            oxs = [oxs; {ox_name, wt, temp, h_specific, exploded_formula}];
        else
            % Prompt the user if they wish to use the ox anyway
            prompt = sprintf('Ox name does not match any reagents in CEA library. Do you wish to use this ox anyway? (y/n): ');
            user_input = input(prompt, 's');

            if strcmpi(user_input, 'y')
                % Input weight fraction
                wt = input('Enter ox weight fraction (note that this value goes from 0 to 1,\nand the sum of all weight fractions if one or more oxs are selected must be equal to 1): ');

                % Redo weight fraction input if invalid input
                while ~(isnumeric(wt) && wt > 0 && wt <= 1)
                    clc;
                    fprintf('Invalid weight fraction value.\n');
                    wt = input('Enter the ox weight fraction: '); 
                end

                % Input reagent temp
                temp = input('Enter the ox temperature in K: ');

                % Redo temp input if invalid input
                while ~(isnumeric(temp) && temp > 0)
                    clc;
                    fprintf('Invalid temperature value.\n');
                    temp = input('Enter the ox temperature in K: '); 
                end

                % Input specific enthalpy in Kj/mol
                h_specific = input('Enter the specific enthalpy in Kj/mol: ');

                % Input exploded chemical formula
                exploded_formula = input('Enter the exploded chemical formula:\nExample: if species is CO2, input should be C 1 O 2\n', 's');

                % Add ox and properties to matrix
                oxs = [oxs; {ox_name, wt, temp, h_specific, exploded_formula}];
            else
                % If the user chooses not to use the ox, prompt again for ox name
                continue;
            end
        end

        % Option to include more oxs
        more = input('Do you wish to add another ox? (y/n): ', 's');
        if ~strcmpi(more, 'y')
            break;
        end
        clc
    end

    % Calculate the sum of weight fractions
    sum_weight_fractions = sum([oxs{:, 2}]);

    % Check if the sum is not equal to 1
    if sum_weight_fractions ~= 1
        fprintf('Error: Sum of weight fractions is not 1. Please input the oxs again.\n');
    else
        flag = 0; % Redo ox inputs if sum of weight fractions is not 1
    end

end



end

