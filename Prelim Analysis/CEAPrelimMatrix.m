function performance_matrix = CEAPrelimMatrix(filename)
% CEAMatrix takes a CEA output file name as an input, and 
% creates an array, where for several values of OF ratio (first column),
% the following columns are the corresponding performance parameters
% associated with that OF ratio value: the second column is the ideal
% nozzle expansion ratio (epsillon), third column is c*, fourth is Thrust
% Coefficient, fifth is vacuum Isp and the sixth is Isp

    % Move to CEA directory where .out files are located
    cd('CEA');

    % Attempt to open the file
    fid = fopen(filename, 'r');
    
    % Check if the file could not be opened
    if fid == -1
        error('Error: Could not open the file.');
    end
    
    % Initialize variables to store data from file
    of_values = [];
    performance_params = [];
    performance = [];
    
    % Read each line of the file
    tline = fgetl(fid);
    while ischar(tline)
        % Check if the line indicates the start of a performance calculation
        if contains(tline, 'COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR')
            % Find the line with O/F value
            while ~contains(tline, 'O/F=')
                tline = fgetl(fid);
            end
            % Extract O/F value
            of_value_parts = regexp(tline, 'O/F=(\s*\d+\.?\d*)', 'tokens');
            if ~isempty(of_value_parts)
                of_value_str = of_value_parts{1}{1};
                of_value = str2double(of_value_str);
                of_values = [of_values; of_value];
                
                % Find the line containing PERFORMANCE PARAMETERS
                while ~contains(tline, 'PERFORMANCE PARAMETERS')
                    tline = fgetl(fid);
                end
                
                % Skip empty line after header
                fgetl(fid);
                
                % Read and extract performance parameters until an empty
                % line is encountered (ie all parameters are recorded)
                tline = fgetl(fid);
                % if line is not empty
                while ~isempty(tline)
                    performance_params_parts = strsplit(tline);
                    % Check if there are at least 4 values in the line
                    if numel(performance_params_parts) >= 4
                        % Extract the 4th numerical value (corresponding to value at nozzle exit)
                        num_values = str2double(performance_params_parts);
                        num_values = num_values(~isnan(num_values));
                        if numel(num_values) >= 2
                            % Extract the second numerical value (corresponding to value at nozzle exit)
                            performance_param = num_values(2);

                            % Save the parameter value in the array
                            performance_params = [performance_params; performance_param];
                        else
                            % Print an error message if there are not enough numerical values
                            fprintf('Error: Unable to compute the performance parameter.\n');
                        end
                    else
                        % If there are not enough values, it indicates an error during CEA analysis
                        % Print an error message
                        fprintf('Error: Unable to compute the performance parameter.\n');
                    end
                    % Read the next line
                    tline = fgetl(fid);
                end
                
                % Records the performance values for a specific OF value in
                % a matrix, and resets performance_params to record
                % parameters for the next OF value
                performance = [performance; performance_params'];
                performance_params = [];
            end
        end
        tline = fgetl(fid);
    end
    
    % Close the file
    fclose(fid);
    
    % Return to MATLAB code directory
    cd('..');
    
    % Combine O/F values and performance parameters into a array, for easy
    % plotting
    performance_matrix = [of_values performance];
end