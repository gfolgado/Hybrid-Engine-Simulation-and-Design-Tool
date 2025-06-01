function data = CEAMatrixv2(filename)
% CEAMatrixv2 takes a CEA output file name as an input, and 
% creates an array, where for several values of OF ratio (first column),
% the following columns are the corresponding performance parameters
% associated with that OF ratio value: 2 column is Pc; 3 column is Pt; 4
% column is P_exit_sub; 5 column is P_exit_sup; 6 column is Tc; 7 column 
% is Tt; 8 columns is rho_c; 9 column is MW_c; 10 column is cp_c; 11 column
% is cp_t; 12 column is gamma_c; 13 column is M_e_sub; 14 column is
% M_e_sup; 15 column is cstar_sub; 16 column is cstar_sup; 17 column is
% CF_sub; 18 column is CF_sup; 19 column is Ivac_sub; 20 column is
% Ivac_sup; 21 column is Isp_sub; 22 column is Isp_sup.

% Move to CEA directory where .out files are located
cd('CEA');
filename = strcat(filename{1}, '.out');

% Attempt to open the file
fid = fopen(filename, 'r');

% Check if the file could not be opened
if fid == -1
    error('Error: Could not open the file.');
    return
end

% Initialize variables to store data from file
of_values = [];
MW_cur = [];
T_cur = [];
P_cur = [];
cp_cur = [];
gamma_cur = 0;
M_cur = [];
data_cur = [];
data_values = [];
cstar_values = [];
Cf_values = [];
Isp_values = [];
performance_params = [];
performance = [];

% Read each line of the file
tline = fgetl(fid);
while ischar(tline)
    % Check if the line indicates the start of a performance calculation
    if contains(tline, 'THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION')
        % Find the line with O/F value
        while ~contains(tline, 'O/F=')
            tline = fgetl(fid);
        end
        % Extract current O/F value
        of_value_parts = regexp(tline, 'O/F=(\s*\d+\.?\d*)', 'tokens');

        if ~isempty(of_value_parts)
            of_value_str = of_value_parts{1}{1};
            of_value = str2double(of_value_str);
            of_values = [of_values; of_value];
        end

        % Find the line containing pressures, and save the values
        while ~contains(tline, 'P, BAR')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        P_strings = line(4:end);
        P_cur = str2double(P_strings);
        % If supersonic conditions are not computed, set value for Pe_sub as
        % -1
        if length(P_cur) == 3
            P_cur(4) = -1;
        end

        % Find line containing temperature, save value for chamber
        while ~contains(tline, 'T, K')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        T_cur = str2double(line(4:5));
        
        % Find line containing density
        while ~contains(tline, 'RHO, KG/CU M')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        % Check if the second last character is + or -
        if line{5}(end-1) == '+' || line{5}(end-1) == '-'
            rho_exponent = str2double(line{5}(end-1:end));
            rho_base = str2double(line{5}(1:end-2));
        else
            rho_exponent = 0;
            rho_base = str2double(line{5});
        end
        rho_cur = rho_base * 10^rho_exponent;       
        % Find line containing Molecular weight
        while ~contains(tline, 'M, (1/n)')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        MW_cur = str2double(line(4));
        MW_cur = MW_cur/1000;
       
        % Find line containing cp value for CC and throat
        while ~contains(tline, 'Cp, KJ/(KG)(K)')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        cp_cur = str2double(line(4:5)) * 10^3;

        % Find line containing gamma for CC
        while ~contains(tline, 'GAMMAs')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        gamma_cur = str2double(line(3));

        % Find line containing Mach numbers for both exit conditions
        while ~contains(tline, 'MACH NUMBER')
            tline = fgetl(fid);
        end
        line = strsplit(tline);
        M_cur = str2double(line(6:end));
        % If M_cur length is 1, there is no value for Me_sup, so it is set
        % as -1
        if length(M_cur) == 1
            M_cur(2) = -1;
        end
        
        data_cur = [P_cur, T_cur, rho_cur, MW_cur, cp_cur, gamma_cur, M_cur];
        
        % Find line containing Performance Parameters
        while ~contains(tline, 'PERFORMANCE PARAMETERS')
            tline = fgetl(fid);
        end

        % Skip empty line after header and Ae/At values that are inputs
        fgetl(fid);
        fgetl(fid);

        % Read and extract performance parameters until an empty
        % line is encountered (ie all parameters are recorded)
        tline = fgetl(fid);
        % if line is not empty
        while ~isempty(tline)
            performance_params_parts = strsplit(tline);
            % Check if there are at least 4 values in the line (if not this
            % implies that CEA could not perform calculations)
            if numel(performance_params_parts) >= 4
                num_values = str2double(performance_params_parts);
                num_values = num_values(~isnan(num_values));
                % Check length of num_values, if it's less than 2, this
                % means that calculations for subsonic case could not be
                % performed, so the value for the performance parameter in
                % the subsonic case is saved as -1
                if length(num_values) == 2
                    performance_param(2) = -1;
                    performance_param(1) = num_values(2);
                else
                    performance_param = num_values(end-1:end);
                end
                % Save the parameter value in the array
                performance_params = [performance_params, performance_param];

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
        % parameters for the next OF value; does the same for other thermo
        % data gathered
        
        data_values = [data_values; data_cur];
        data_cur  = [];
        T_cur = [];
        P_cur = [];
        rho_cur = [];
        MW_cur = [];
        cp_cur = [];
        M_cur = [];
        
        performance = [performance; performance_params];
        performance_params = [];
    end
    tline = fgetl(fid);
end
    

% Close the file
fclose(fid);

% Return to MATLAB code directory
cd('..');

% Combine O/F values and performance parameters into a array, for easy
% plotting
data = [of_values data_values performance];
end