%% Hybrid Engine Simulation Tool v4

%% Intro/Requirements

clc
clear all

% Master's Degree Thesis in Aerospace Engineering
% Script developed by Guilherme das Neves Pereira Navarro Folgado, 92685

% This script is the final amalgamation of most of the developed
% functions/scripts for the oxidiser tank, injector and combustion chamber
% dynamics.
% The code uses a CEA automation tool to find combustion properties, and
% manipulated N2O saturation data (TankData struct) alongside the CoolProp
% software to find the properties of the oxidiser along several points of
% the feed system. As such, the folder labeled CEA should not be deleted or
% altered, and a distribution of Python is required in order to run
% CoolProp (recommended Python v3.6, which despite not being supported, is
% still avaliable online and is compatible with CoolProp).

% Most if not all of the required inputs from the engine design can be
% estimated initially using the developed fuel dimensioning script entitled
% "FuelDimensioning_NozzleDesign.m" and by providing some desired
% performance parameters (thrust, expected burntime, initial inner core
% radius).

%% Select Loading Simulation Data or Running Simulation

while true
    user_choice = input('Do you wish to load previous simulation data or run a new simulation? (type "load" or "run"): ', 's');
    if strcmpi(user_choice, 'load') || strcmpi(user_choice, 'run')
        break;
    else
        clc
        disp('Invalid choice. Please type "load" or "run".');
    end
end

%% Loading Data, Tables and Variables

load('N2OSaturatedPropertiesv2.mat');
N2OSat.metadata = metadata;
N2OSat.dataset = str2double(N2OData);
load('TankConditions4.mat');
load('TankProperties4.mat');
TankData.TankConditions = TankConditions4;
TankData.TankProperties = TankProperties4;

if strcmpi(user_choice, 'load')
    % Check for the existence of the "Results" folder
    if ~isfolder('Results')
        error('Could not locate the Results folder.');
    end
    
    % List available simulation data
    results_folders = dir('Results');
    results_folders = results_folders([results_folders.isdir] & ~startsWith({results_folders.name}, '.'));
    
    if isempty(results_folders)
        error('No simulation data found in the Results folder.');
    end
    
    disp('Available simulation data:');
    for i = 1:length(results_folders)
        fprintf('%d: %s\n', i, results_folders(i).name);
    end
    
    % Prompt user to select a folder
    while true
        selected_folder = input('Which simulation data shall be loaded? (type the index or folder name): ', 's');
        
        if ~isempty(str2double(selected_folder)) && ...
                str2double(selected_folder) > 0 && str2double(selected_folder) <= length(results_folders)
            selected_folder = results_folders(str2double(selected_folder)).name;
        end
        
        if isfolder(fullfile('Results', selected_folder))
            break;
        else
            disp('Invalid selection. Please type a valid index or folder name.');
        end
    end
    
    % Load data from the selected folder
    data_file = fullfile('Results', selected_folder, 'simulation_data.mat');
    if isfile(data_file)
        load(data_file, 'results', 'time');
        disp(['Loaded simulation data from folder: ', selected_folder]);
    else
        error('Simulation data filea not found in the selected folder.');
    end
    
else
    
    %% User Input for Model Selection
    
    models = {'spi', 'hem', 'dyer', 'spc', 'fml'}; % Injector mass flow models
    
    while true
        disp('Available injector mass flow models:');
        for j = 1:length(models)
            fprintf('%d: %s\n', j, models{j});
        end
        fprintf('%d: Run all models\n', length(models) + 1);
        
        choice = input('Select a model to run (index): ', 's');
        choice_num = str2double(choice);
        
        clc
        
        if ~isnan(choice_num) && choice_num >= 1 && choice_num <= length(models) + 1
            if choice_num == length(models) + 1
                selected_models = models;
            else
                selected_models = models(choice_num);
            end
            break;
        else
            clc
            disp('Invalid input. Please select a valid model index.');
        end
    end
    
    results = struct();     % Struct to save results for different mass flow models
    
    %% Initializing Variables

    tic
    
    for model_idx = 1:length(selected_models)
        model = selected_models{model_idx};
        
        clearvars -except results model models selected_models model_idx TankData N2OSat;
        
        % Tank Simulation Variables
        
        MW = 0.044013;          % Molecular weight of N2O [kg/mol]
        V_tank = 0.03;          % Tank Volume [m^3]
        m_tank(1) = 19.96;      % Initial oxidiser mass in tank [kg]
        P_tank(1) = 6;          % Initial Tank Pressure [MPa]
        delta_p = 0.5;          % Feed System Pressure Losses [MPa]
        Cd_inj = 0.63;          % Injector Discharge Coefficient
        Cd_vent = 0.7;          % Vent Discharge Coefficient
        A_inj = 3.1415926e-5;   % Injector Total Orifice Area [m^2]
        A_vent = 5e-7;          % Vent Total Orifice Area [m^2]
        tstop = 60;             % Max Simulation TIme [s]
        tstep = 0.25;            % Simulation Time Step [s]
        lro = NaN;              % Liquid Run out Time [s]
        endtime = NaN;          % Simulation Finish Time [s]
        
        % Constants Struct for Mass Flow Function input
        const.A_inj = A_inj;
        const.Cd_inj = Cd_inj;
        const.delta_p = delta_p;
        const.N2OSat = N2OSat;
        const.MW = MW;
        
        % Combustion Simulation Variables
        
        rho_f = 900;                % Fuel grain density [kg/m^3]
        a = 0.155 * 10^-3;          % Regression rate coefficient
        n = 0.5;                    % Regression rate exponent
        L_g = 0.3;                  % Fuel grain length [m]
        zeta_d = 1.07;              % Nozzle discharge correction factor
        zeta_cstar = 0.95;           % Characteristic velocity correction factor
        zeta_CF = 0.95;              % Thrust Coefficient correction factor
        r_chamber = 0.156/2;        % Combustion Chamber Radius [m]
        A_ratio_nozzle = 5.99;      % Nozzle area ratio
        r_t = 0.0298/2;             % Nozzle Throat Radius [m]
        r_e = 0.0731/2;             % Nozzle Exit Radius [m]
        T_amb = 298.15;             % Ambient Temperature [K]
        p_amb = 1.5 * 10^5;      % Ambient Pressure [Pa]
        r_p_init = 0.02;            % Fuel Grain Port Initial Radius [m]
        iteration_tolerance = .01;  % 1% tolerance for iterative solving
        R = 8.314462;               % Universal Gas Constant
        
        % Combustion Chamber Model Variables Initialization
        
        A_t = pi() * (r_t)^2; % Nozzle throat area [m^2]
        A_ratio_nozzle_eps = iteration_tolerance * A_ratio_nozzle; % Acceptable error of nozzle area ratio
        fuels = {'C32H66(a)', 0.95, 298.15, 0, ''; 'stacid', 0.05, 298.15, -912, 'C 18 H 36 O 2'}; % List of fuels for CEA
        oxs = {'N2O(L),298.15K', 1, 298.15, 0, ''}; % List of oxidizers for CEA
        m_tot = 0;          % Total propellant mass in CC [kg]
        r_p = r_p_init;     % Initial fuel grain port radius [m]
        p_c = p_amb;        % Chamber pressure initialized as ambient pressure [Pa]
        T_c = T_amb;        % Chamber temperature initialized as ambient temp [K]
        Tstag_c = T_amb;    % Chamber stagnation temperature initialized as ambient temp [K]
        M_air = 0.02897;    % Air molar mass at 298K [kg/mol]
        R_c = R / M_air;    % Specific gas constant, initialized for air since no combustion yet [J/kg*K]
        k_c = 1.4;          % Specific heat ratio; initialized for air at 298K
        Ma = 3;             % Exit Mach number
        burntime = tstop;   % Burn duration, initialized as max simulation time [s]
        const.Pc = p_c;
        
        %% Initial Conditions Calculations
        
        % Initial Fluid Density in tank [kg/m^3]
        rho_tank(1) = m_tank(1)/V_tank;
        
        % Find intial conditions in tank
        data = InitConditions(TankData, P_tank(1), rho_tank(1), 'P');
        
        s_tank(1) = data(1);    % Initial Tank Specific Entropy [J/kg.K]
        x(1) = data(2);         % Initial Tank Fluid Quality
        P_tank(1) = data(4);    % Initial Tank Pressure [MPa]
        T_tank(1) = data(3);    % Initial Tank Temperature [K]
        h_tank(1) = data(5);    % Initial Tank Specific Enthalpy [kJ/kg]
        
        if x(1)<0 || x(1)>1
            error('Computed initial tank fluid quality is impossible. Check the script inputs');
        end
        
        data = N2OSatInterpolator (N2OSat, T_tank(1), 'T', {'rho_v','rho_l','s_v','s_l'});
        
        rho_v(1) = data(1);         % Vapour density [kg/m^3]
        rho_l(1) = data(2);         % Liquid density [kg/m^3]
        s_v(1) = data(3)/MW;        % Vapour specific entropy [J/kg.K]
        s_l(1) = data(4)/MW;        % Liquid specific entropy [J/kg.K]
        
        % Initial Oxidiser Mass in tank [kg]
        m_l(1) = (1 - x(1)) * m_tank(1);
        m_v(1) = x(1) * m_tank(1);
        
        if m_l(1) < 0 || m_v(1) < 0
            error('Values for liquid/vapour phase mass in tank are negative. Double check inputs.');
        end
        
        % Initial Tank Total Entropy [J/K]
        S_tank(1) = s_tank(1) * m_tank;
        
        %% Simulation Loop
        
        i_f = tstop/tstep;      % Final Simulation timestep for input tstop value
        
        for i = 1:i_f
            const.i = i;
            %% Tank Model
            
            % Simulation Stop Conditions (nitrous in tank depleted, pressure
            % differential between tank and chamber insufficient to sustain flow
            
            fprintf('Timestep %f\n', i);
            
            if m_tank(i) < 0.05
                fprintf('No more nitrous left in tank\n');
                endtime = (i-1) * tstep;
                break
            end
            
            if P_tank(i)*10^6 - delta_p*10^6 - p_c(i) < 1e4
                fprintf('Pressure difference between tank and chamber is no longer sufficient to sustain oxidiser flow. Ending simulation.\n');
                endtime = (i - 1) * tstep;
                break
            end
            
            % Liquid Runout Detection
            
            if x(i) >= 1 && isnan(lro)
                lro = (i-2) * tstep;
                fprintf('Liquid Runout at %.2f seconds\n', lro);
            end
            
            % Update const struct with new values
            const.x = x(i);
            const.rho_l = rho_l(i);
            
            % Oxidiser Mass Flow Rate into combustion chamber [kg/s]
            if isnan(lro)
                mdot_ox(i) = MassFlow(TankData, P_tank(i), const, s_tank(i), rho_tank(i), model);
            else
                gamma = Cp_v / Cv_v;        % Ratio of Specific Heats
                mdot_ox(i) = A_inj  * Cd_inj * sqrt(gamma * P_tank(i) * 10^6 * rho_v(i) * (2 / (gamma + 1))^((gamma + 1)/(gamma - 1)));
            end
            
            % Oxidiser Mass Flow Rate through vent [kg/s]
            mdot_vent(i) = A_vent * Cd_vent * sqrt(2 * rho_v(i) * (P_tank(i) * 10^6 - p_amb));
            
            if i == 1  % Take average mass flow over last two time periods to attenuate numerical instability
                mdot_ox(i) = 1/2*mdot_ox;
                mdot_vent(i) = 1/2 * mdot_vent;
            else
                mdot_ox(i) = 1/2 * mdot_ox(i) + 1/2 * mdot_ox(i-1);
                mdot_vent(i) = 1/2 * mdot_vent(i) + 1/2 * mdot_vent(i-1);
            end
            
            % Total mass leaving tank during timestep [kg]
            dm_ox = mdot_ox(i) * tstep;
            dm_vent = mdot_vent(i) * tstep;
            
            %% Combustion Model
            
            A_burn(i) = pi() * r_p(i)^2; % Fuel grain burning area [m^2]
            G(i) = mdot_ox(i) / A_burn(i); % Oxidiser Mass Flux [kg/s*m^2]
            dr_p(i) = a * G(i)^n * tstep; % Regression rate [m/s]
            mdot_f(i) = 2 * pi() * r_p(i) * L_g * rho_f * (dr_p(i) / tstep); %Fuel mass flow rate [kg/s]
            
            % Iteratively find actual fuel mass flow value
            k = 0;
            mdot_f_temp = 0;
            while abs(mdot_f_temp - mdot_f(i)) > 0.01*mdot_f(i) && k < 100
                mdot_f_temp = mdot_f(i);
                mdot_cc(i) = mdot_f(i) + mdot_ox(i); % Total gaseous propellant mass flow into CC [kg/s]
                G(i) = (mdot_ox(i) + mdot_cc(i)) / (2 * A_burn(i)); % Average mass flux in fuel grain port [kg/s*m^2]
                dr_p(i) = a * G(i)^n * tstep; % Regression rate for updated mass flux [m/s]
                mdot_f(i) = 2 * pi() * r_p(i) * L_g * rho_f * (dr_p(i) / tstep); % Fuel mass flow rate [kg/s]
                k = k + 1;
            end
            
            OF(i) = mdot_ox(i) ./ mdot_f(i); % Oxidizer to fuel ratio
            % Steady state choked flow expression for pressure
            pstag_c(i) = mdot_cc(i) / (zeta_d * A_t) * sqrt(Tstag_c(i) * ...
                R_c(i) / k_c(i) * ((k_c(i) +1) / 2)^((k_c(i) + 1) / (k_c(i) - 1)));
            
            % Chamber flow velocity calculation (first time step is ignored bc its too cold)
            if i > 1
                v_c(i) = mdot_cc(i) / (rho_c(i) * A_burn(i));
                v_c(i) = 0.5 * v_c(i) + 0.5 * v_c(i-1); % Do average of velocity with previous time step to attenuate numerical instabilities
                T_c(i) = Tstag_c(i) - v_c(i)^2 / (2 * cp_c(i));
            else
                T_c(i) = Tstag_c(i);
            end
            
            p_c(i) = pstag_c(i) * ((T_c(i) / Tstag_c(i)) .^ (k_c(i) / (k_c(i) - 1)));   % Compute chamber pressure with velocity corection
            
            if i >= 2
                p_c(i) = 1/2 * (p_c(i) + p_c(i-1));
            end
            
            [dataProp, dataNozzle] = CEAProp(p_c(i)/10^5, OF(i), A_ratio_nozzle, fuels, oxs);  % Get CEA combustion data
            dataProp = num2cell(dataProp);
            [Tstag_c(i+1), rho_c(i+1), cp_c(i+1), k_c(i+1), M_c(i+1)] = ...
                dataProp{:};
            
            Tstag_c(i+1) = Tstag_c(i+1) * zeta_cstar^2; % Correct stagnation temperature
            R_c(i+1) = R / M_c(i+1);  % Update specific gas constant
            
            Ain.k = k_c(i);
            Ain.A = A_ratio_nozzle;
            if abs(Aerror(Ma(i), Ain)) > A_ratio_nozzle_eps
                Ma(i) = secant(@Aerror, Ma(i), Ain);
            end
            
            p_e(i) = pstag_c(i) ./ (1 + (k_c(i) - 1)./2 .* Ma(i).^2).^(k_c(i) / (k_c(i) -1));
            
            T_e(i) = Tstag_c(i) ./ (1 + (k_c(i) - 1) ./ 2 .* Ma(i).^2); % Flow temperature at nozzle exit [K]
            v_e(i) = Ma(i) .* sqrt(k_c(i) .* R_c(i) .* T_e(i)); % Fluid velocity at nozzle exit [m/s]
            A_ratio_eff = A_ratio_nozzle;

            F(i) = zeta_CF * (mdot_cc(i) .* v_e(i) + (p_e(i) - p_amb) .* A_t .* A_ratio_eff); % Thrust generated

            mdot_nozzle(i) = F(i)/v_e(i);

            %% Update Tank/Combustion Model Variables
            
            % Iterate in time if max sim time has not been reached, and
            % there is still more than 5% of initial oxidiser mass in the
            % tank, and if the fuel grain radius has not exceeded the CC radius
            if (i < tstop/tstep) && (m_tank(i) / m_tank(1) > 0.05) && (r_chamber - r_p(i) > 0)
                m_tot(i+1) = m_tot(i) - mdot_cc(i) * tstep;
                r_p(i+1) = r_p(i) + dr_p(i);
                p_c(i+1) = p_c(i);
                Ma(i+1) = Ma(i);
            else
                if (m_tank(i) / m_tank(1) < 0.05)
                    fprintf('Oxidiser mass remaining inside the tank is lower than 5%% of the starting value. Ending simulation.\n');
                end
                if (r_chamber - r_p(i) < 0)
                    fprintf('Fuel grain has been entirely combusted. Ending simulation.\n');
                end
                endtime = (i - 1) * tstep;
                break
            end
            
            if isnan(lro)
                S_tank(i+1) = S_tank(i) - (dm_ox * s_l(i)) - (dm_vent * s_v(i));    % Total Tank Entropy [J/K]
            else
                S_tank(i+1) = S_tank(i) - ((dm_ox + dm_vent) * s_v(i));
            end
            m_tank(i+1) = m_tank(i) - dm_ox - dm_vent;                          % Total Oxidiser mass [kg]
            s_tank(i+1) = S_tank(i+1) / m_tank(i+1);                            % Tank Specific Entropy [J/kg.K]
            rho_tank(i+1) = m_tank(i+1) / V_tank;                               % Tank Fluid Density [kg/m^3]
            
            % Find next time step tank conditions from manipulated N2O Saturation data
            out = TankPropsInterpolator (TankData, s_tank(i+1)/10^3, rho_tank(i+1), {'x','P','T','h'});
            
            x(i+1) = out(1);
            if ~isnan(lro)
                x(i+1) = 1;
            end
            T_tank(i+1) = out(2);
            P_tank(i+1) = out(3);
            h_tank(i+1) = out(4);
            
            % Find relevant tank fluid properties via interpolation of N2O Saturation dataset from NIST
            out = N2OSatInterpolator (N2OSat, P_tank(i+1), 'P', {'rho_v','rho_l','s_v','s_l','Cv_v','Cp_v'});
            
            rho_v(i+1) = out(1);
            rho_l(i+1) = out(2);
            s_v(i+1) = out(3)/MW;
            s_l(i+1) = out(4)/MW;
            Cv_v = out(5);
            Cp_v = out(6);
            
            % Update Tank vapour and liquid mass from total tank mass and fluid
            % quality
            m_v(i+1) = x(i+1) * m_tank(i+1);
            m_l(i+1) = (1 - x(i+1)) * m_tank(i+1);
            const.Pc = p_c(i);
        end
        
        if isnan(endtime)
            lro = tstop;
            endtime = tstop;
        end
        
        % Store results for current mass flow model
        results.(model).endtime = endtime;
        results.(model).P_tank = P_tank;
        results.(model).T_tank = T_tank;
        results.(model).m_tank = m_tank;
        results.(model).rho_tank = rho_tank;
        results.(model).x = x;
        results.(model).s_tank = s_tank;
        results.(model).S_tank = S_tank;
        results.(model).h_tank = h_tank;
        results.(model).m_l = m_l;
        results.(model).m_v = m_v;
        results.(model).mdot_ox = mdot_ox;
        results.(model).mdot_vent = mdot_vent;
        results.(model).F = F;
        results.(model).OF = OF;
        results.(model).p_c = p_c;
        results.(model).T_c = T_c;
        results.(model).p_e = p_e;
        results.(model).T_e = T_e;
        results.(model).r_p = r_p;
        results.(model).lro = lro;
    end
    
    toc
    
    %% Time Vector Setup
    
    for model = fieldnames(results)'
        model_name = model{1};
        if isfield(results.(model_name), 'endtime')
            model_endtime = results.(model_name).endtime;
            if model_endtime > endtime
                endtime = model_endtime;
            end
        end
    end
    
    time = 0:tstep:endtime; % Create the time vector
    
    %% Saving Simulation Data
    
    % After simulation, prompt to save data
    while true
        save_choice = input('Do you wish to save the current simulation data? (yes or no): ', 's');
        if strcmpi(save_choice, 'yes') || strcmpi(save_choice, 'no')
            break;
        else
            disp('Invalid choice. Please type "yes" or "no".');
        end
    end
    
    if strcmpi(save_choice, 'yes')
        
        % Create the Results folder if it does not exist
        if ~isfolder('Results')
            mkdir('Results');
        end
        
        while true
            save_folder = input('Enter the save folder name: ', 's');
            save_path = fullfile('Results', save_folder);
            
            if isfolder(save_path)
                overwrite = input('Folder already exists. Do you wish to overwrite? (yes or no): ', 's');
                if strcmpi(overwrite, 'yes')
                    break;
                elseif strcmpi(overwrite, 'no')
                    disp('Please enter a new folder name.');
                    continue;
                else
                    disp('Invalid input. Please respond with yes or no.');
                end
            else
                mkdir(save_path);
                break;
            end
        end
        
        % Save the data
        results_folder_path = fullfile('Results', save_folder);
        save(fullfile(results_folder_path, 'simulation_data.mat'), 'results', 'time');
        disp(['Simulation data saved in folder: ', save_folder]);
    else
        disp('Simulation data not saved.');
    end
end

%% Plotting Loop

% Define available variables for plotting
variables = {'P_tank', 'T_tank', 'm_tank', 'rho_tank', 'x', 's_tank', 'S_tank', 'h_tank', 'm_l', 'm_v', 'mdot_ox', 'mdot_fuel', 'mdot_vent', 'F', 'OF', 'p_c', 'T_c', 'p_e', 'T_e', 'r_p'};

while true
    % Prompt user for variable to plot
    disp('Available variables to plot:');
    for i = 1:length(variables)
        fprintf('%d: %s\n', i, variables{i});
    end
    user_input = input('Enter the variable name or index to plot (type "STOP" to exit): ', 's');

    % Check for termination condition
    if strcmpi(user_input, 'STOP')
        disp('Exiting plotting loop.');
        break;
    end

    % Validate user input
    if ismember(user_input, variables)
        var_to_plot = user_input;
    elseif ~isempty(str2double(user_input)) && str2double(user_input) > 0 && str2double(user_input) <= length(variables)
        var_to_plot = variables{str2double(user_input)};
    else
        disp('Invalid input. Please enter a valid variable name or index.');
        continue;
    end

    % Plot the selected variable for all models
    figure;
    hold on;
    legend_entries = {};
    model_colors = lines(length(fieldnames(results))); % Generate distinct colors
    model_idx = 1;
    plot_handles = [];
    
    for model = fieldnames(results)'
        model_name = model{1};
        if isfield(results.(model_name), var_to_plot)
            plot_data = results.(model_name).(var_to_plot);
            color = model_colors(model_idx, :);
            h = plot(time(1:length(plot_data)), plot_data, 'Color', color, 'DisplayName', model_name, 'LineStyle', '-');
            plot_handles = [plot_handles, h];
            legend_entries{end+1} = model_name;
            
            % Compute and plot average value
            avg_value = mean(plot_data);
            plot([time(1), time(end)], [avg_value avg_value], 'Color', color, 'LineStyle', '--', 'HandleVisibility', 'off');
            
            % Plot liquid runout line if valid
            lro_time = results.(model_name).lro;
            if ~isnan(lro_time)
                yl = ylim;
                plot([lro_time, lro_time], yl, '--', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
                % Add vertical text label
                text(lro_time, yl(2), sprintf('%s runout at %.1f s', model_name, lro_time), ...
                    'Rotation', 90, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
                    'Color', color, 'FontSize', 8);
            end
            
            % If the variable is 'F', calculate and display total impulse
            if strcmp(var_to_plot, 'F')
                total_impulse = trapz(time(1:length(plot_data)), plot_data);
                fprintf('Total Impulse: %.2f N.s\n', total_impulse);
            end
            
            model_idx = model_idx + 1;
        end
    end
    hold off;

    % Customize plot
    title(['Plot of ', var_to_plot, ' with time']);
    xlabel('Time [s]');
    ylabel(var_to_plot);
    legend(plot_handles, legend_entries, 'Location', 'best');
    grid on;
end

%% Subfunctions Definition

function [out] = N2OSatInterpolator (N2OSat ,target, in_type, out_type)
% SatInterpolator computes the saturation temperature N2O properties for
% the selected input values that may not be a defined data point in the
% data table
% [out] - numeric array containing the computed saturation thermodynamic
% properties
% N2OSat - struct containing the data points for saturated N2O and the
% associated metadata
% target - numeric input value for which the thermodynamic fluid properties
% are to be found
% in_type - string defining the input thermodynamic property (pressure,
% temp, entropy, etc.)
% out_typ - array of strings defining the desired thermodynamic outputs to
% be calculated for the speficified input value

    out = zeros(1,length(out_type));

    % Identify column index for selected input type
    in_col = find(strcmp(N2OSat.metadata(2, :), in_type), 1);
    if isempty(in_col)
        error('Input type "%s" not found in metadata.', in_type);
    end

    % Identify column indices for selected output types
    out_cols = zeros(1, length(out_type));
    for i = 1:length(out_type)
        col = find(strcmp(N2OSat.metadata(2, :), out_type{i}));
        if isempty(col)
            error('Output type "%s" not found in metadata.', out_type{i});
        end
        out_cols(i) = col;
    end

    % Find surrounding points
    idx_below = find(N2OSat.dataset(:, in_col) <= target, 1, 'last');
    idx_above = find(N2OSat.dataset(:, in_col) >= target, 1, 'first');

    if isempty(idx_below) || isempty(idx_above)
        error('Target value %.2f is out of the data range.', target);
    end

    if idx_below == idx_above
        % Exact match, no interpolation needed
        out(:) = N2OSat.dataset(idx_below, out_cols);
    else
        % Linear interpolation
        x1 = N2OSat.dataset(idx_below, in_col);
        x2 = N2OSat.dataset(idx_above, in_col);
        y1 = N2OSat.dataset(idx_below, :);
        y2 = N2OSat.dataset(idx_above, :);

        interpolated_row = y1 + (y2 - y1) * (target - x1) / (x2 - x1);
        out(:) = interpolated_row(out_cols);
    end
end


function [out] = TankPropsInterpolator (TankData, s_target, rho_target, out_type)
% TankPropsInterpolator is able compute the N2O fluid properties in the
% oxidiser tank, by interpolating a known data set for a target value
% of specific entropy and density

% TankData - struct containing a 3D matrix named TankProperties, which
% holds the dataset to be interpolated, where every 3rd dimension array
% is a different property (in order, the fluid quality, the tank
% temperature, the tank pressure and the specific enthalpy). This
% struct also holds the TankConditions which is a struct in and of
% itself, and it contains two arrays, one named S and one named rho.
% The S array is a column representing the specific entropy values for
% which the TankProperties were computed, and the rho array is a row
% array containing the values for rho for which the results were
% computed.
% s_target - numerical value for the specific entropy for which the
% interpolation should be performed
% rho_target - numerical value for the tank density for which the
% interpolation should be performed
% out_type - array of strings containing the desired outputs of the
% interpolation ('x' for quality, 'P' for pressure, 'T' for temperature
% and 'h' for specific enthalpy)
% out - numeric array containing the desired outputs for the target s
% and rho values

    % Initialize output array
    out = zeros(1, length(out_type));

    % Extract results and associated values of s and rho
    S = TankData.TankConditions.S; % Column array for specific entropy
    rho = TankData.TankConditions.rho; % Row array for density
    TankProperties = TankData.TankProperties; % 3D matrix for properties

    % Find indices surrounding s_target, compute fraction for interpolating
    if s_target < S(1) || s_target > S(end)
        error('s input value %.3f is outside the data range.', s_target);
    end
    s_low_idx = find(S <= s_target, 1, 'last');
    s_high_idx = find(S >= s_target, 1, 'first');
    if s_low_idx == s_high_idx
        s_frac = 0;
    else
        s_frac = (s_target - S(s_low_idx)) / (S(s_high_idx) - S(s_low_idx));
    end

    % Find indices surrounding rho_target, compute fraction for
    % interpolating
    if rho_target < rho(1) || rho_target > rho(end)
        error('rho_target value %.3f is outside the interpolation range.', rho_target);
    end
    rho_low_idx = find(rho <= rho_target, 1, 'last');
    rho_high_idx = find(rho >= rho_target, 1, 'first');
    if rho_low_idx == rho_high_idx
        rho_frac = 0;
    else
        rho_frac = (rho_target - rho(rho_low_idx)) / (rho(rho_high_idx) - rho(rho_low_idx));
    end

    % Perform interpolation and save the requested outputs in the out array
    for i = 1:length(out_type)
        % Find 3rd dimension index for extracting desired datapoints
        switch out_type{i}
            case 'x', prop_idx = 1; % Fluid quality
            case 'T', prop_idx = 2; % Temperature
            case 'P', prop_idx = 3; % Pressure
            case 'h', prop_idx = 4; % Specific enthalpy
            otherwise
                error('Invalid property in out_type: %s', out_type(i));
        end

        % Get values at the correct indices
        Q11 = TankProperties(s_low_idx, rho_low_idx, prop_idx);
        Q12 = TankProperties(s_low_idx, rho_high_idx, prop_idx);
        Q21 = TankProperties(s_high_idx, rho_low_idx, prop_idx);
        Q22 = TankProperties(s_high_idx, rho_high_idx, prop_idx);

        % Perform bilinear interpolation
        Q_interp = (1 - s_frac) * (1 - rho_frac) * Q11 + ...
            (1 - s_frac) * rho_frac * Q12 + ...
            s_frac * (1 - rho_frac) * Q21 + ...
            s_frac * rho_frac * Q22;

        % Save interpolated value to output array
        out(i) = Q_interp;
    end
end


function [out] = InitConditions(TankData, PT_in, rho_tank, in_type)
% InitConditions provides the initial tank conditions (namely the tank
% fluid quality and specific entropy) from user input values (either
% pressure or temperature), and from the computed initial tank fluid
% density (from loaded prop mass and tank volume)

% TankData - struct containing a 3D matrix named TankProperties, which
% holds the dataset to be interpolated, where every 3rd dimension array
% is a different property (in order, the fluid quality, the tank
% temperature, the tank pressure and the specific enthalpy). This
% struct also holds the TankConditions which is a struct in and of
% itself, and it contains two arrays, one named S and one named rho.
% The S array is a column representing the specific entropy values for
% which the TankProperties were computed, and the rho array is a row
% array containing the values for rho for which the results were
% computed.
% PT_in - numerical input for either the tank initial temperature [K]
% or pressure [MPa]
% rho_tank - numerical input for the initial tank fluid density,
% computed from the inputs for loaded mass and tank volume [kg/m^3]
% in_type - string to differentiate the type of the input for PT_in.
% Can either be 'P' for pressure input or 'T' for temperature input.

% out - numerical array containing the relevant initial tank conditions
% which are in order: s_tank [J/kg.K], fluid quality, pressure [MPa]
% and temperature [K]

    rho = TankData.TankConditions.rho;
    S = TankData.TankConditions.S;
    TankProperties = TankData.TankProperties;

    switch in_type
        case 'P'
            idx_inp = 2;
        case 'T'
            idx_inp = 3;
        otherwise
            error('Invalid input type in function InitConditions. Must be either P or T');
    end

    if rho_tank < rho(1) || rho_tank > rho(end)
        error('Density input value %.3f is outside the data range.', rho_tank);
    end

    rho_low_idx = find(rho <= rho_tank, 1, 'last');
    rho_high_idx = find(rho >= rho_tank, 1, 'first');

    x1 = rho(rho_low_idx);
    x2 = rho(rho_high_idx);
    y1 = TankProperties(:,rho_low_idx, idx_inp);
    y2 = TankProperties(:,rho_high_idx, idx_inp);

    interpolated_array = y1 + (y2 - y1) * (rho_tank - x1) / (x2 - x1);

    PT_low_idx = find(interpolated_array <= PT_in, 1, 'last');
    PT_high_idx = find(interpolated_array >= PT_in, 1, 'first');

    if isempty(PT_high_idx) ||  isempty(PT_low_idx)
        error('Initial Tank Conditions could not be computed. This likely means N2O cannot be saturated at the provided conditions. Please alter the inputs.');
    end

    x1 = interpolated_array(PT_low_idx);
    x2 = interpolated_array(PT_high_idx);
    y1 = S(PT_low_idx);
    y2 = S(PT_high_idx);

    s_tank = y1 + (y2 - y1) * (PT_in - x1) / (x2 - x1);

    out(1) = s_tank * 10^3;
    holder = TankPropsInterpolator (TankData, s_tank, rho_tank, {'x','P','T', 'h'});

    out = cat(2,out,holder);
end

function [out] = MassFlow (TankData, P_tank, const, s_tank, rho_tank, model)
% MassFlow function - computes the mass flow from the oxidiser tank into
% the combustion chamber via the injector assembly. Can compute this value
% for several different mass flow models found in literature. This specific
% funcion requires the use of CoolProp software to find derivatives of
% fluid properties (exclusive for SPC and FML models).

% TankData - struct containing data manipulated from N2O saturated data
% found in NIST Webbook. Allows quick computation of tank conditions such
% as pressure, temperature, specific enthalpy and fluid quality via the
% TankDataInterpolator function. Used to find relevant data to compute mass
% flow;
% P_tank - oxidiser tank pressure in MPa;
% const - struct containing several required constant values for computing
% the mass flow, such as the injector oriffice area, discharge coefficient,
% injector pressure drop, N2O molecular weight, etc;
% s_tank - tank specific entropy in J/K.kg;
% rho_tank - tank fluid density in kg/m^3;
% model - string indicating which of the desired mass flow models should be
% used to compute the output. Can be either SPI, SPC, HEM, Dyer or FML
% model.

% out - numerical output for mass flow computed using the selected mass
% flow model. Result in kg/s.

    rho = TankData.TankConditions.rho;
    S = TankData.TankConditions.S * 10^3;
    TankProperties = TankData.TankProperties;

    A_inj = const.A_inj;        % Injector Oriffice Area [m^2]
    Cd_inj = const.Cd_inj;      % Injector Discharge Coefficient
    delta_p = const.delta_p;    % Injector Pressure Drop [MPa]
    rho_l = const.rho_l;        % Liquid Phase Nitrous Density [kg/m^3]
    Pc = const.Pc;              % Combustion Chamber Pressure [MPa]
    N2OSat = const.N2OSat;      % Saturated N2O Data
    MW = const.MW;              % N2O Molecular Weight
    R = 8.314462;               % Gas Constant [J/mol.K]
    i = const.i;

    if strcmp(model, 'spi')
        P1 = P_tank - delta_p;
        
        out = A_inj * Cd_inj * sqrt(2 * rho_l * (P1 * 10^6 - Pc));

    elseif strcmp(model, 'hem_calc')
        
        h1 = TankPropsInterpolator(TankData, s_tank/10^3, rho_tank, {'h'}) * 10^3;
        
        aux_v = num2cell(PropsSI({'S', 'H'}, 'P', Pc, 'Q', 1, 'N2O'));
        aux_l = num2cell(PropsSI({'S', 'H'}, 'P', Pc, 'Q', 0, 'N2O'));

        [s_v, h_v] = aux_v{:};
        [s_l, h_l] = aux_l{:};
        
        x2 = (s_tank - s_l)/(s_v - s_l);
                
        h2 = x2 * h_v + (1-x2) * h_l;
        
        rho2 = PropsSI('D', 'P', Pc, 'H', h2, 'N2O');
        
        out = A_inj * Cd_inj * rho2 * sqrt(2 * (h1 - h2));
        
    elseif strcmp(model, 'hemc')

        mdot = 0;
        for P2 = 0.1:0.05:P_tank
            const_temp = const;
            const_temp.Pc = P2 * 10^6;
            aux = MassFlow(TankData, P_tank, const_temp, s_tank, rho_tank, 'hem_calc');
            if aux > mdot
                mdot = aux;
                P_crit = P2;
            else
                break 
            end
        end
        
        out.Pcrit = P_crit;
        out.mdotc = mdot;
        
    elseif strcmp(model, 'hem')
        
        P1 = P_tank - delta_p; 
        
        crit = MassFlow(TankData, P1, const, s_tank, rho_tank, 'hemc');
        P_crit = crit.Pcrit * 10^6;
        
        if Pc < P_crit
            out = crit.mdotc;
        else
            out = MassFlow(TankData, P1, const, s_tank, rho_tank, 'hem_calc');
        end

    elseif strcmp(model, 'dyer')
        
        % Upstream Pressure [MPa]
        P1 = P_tank - delta_p;
        % Downstream pressure [MPa]
        P2 = Pc / 10^6;
        % Non-Equilibrium Coefficient
        k = sqrt((P1*10^6 - P2*10^6)/(P1*10^6 - P2*10^6));

        % Mass flows for SPI and HEM models
        mdot_SPI = MassFlow (TankData, P_tank, const, s_tank, rho_tank, 'spi');
        mdot_HEM = MassFlow (TankData, P_tank, const, s_tank, rho_tank, 'hem');

        % Dyer Mass Flow Model equation
        out = (mdot_HEM + k * mdot_SPI)/(1 + k);

    elseif strcmp(model, 'spc')
        % Upstream Pressure [MPa]
        P1 = P_tank - delta_p;
        % Downstream pressure [MPa]
        P2 = const.Pc/10^6;

        % Find relevant tank data
        data = N2OSatInterpolator (N2OSat, P_tank, 'P', {'rho_l','T','Cp_l','Cv_l'});

        rho_l1 = data(1);   % Liquid phase density [kg/m^3]
        T_t = data(2);      % Tank Temperature [K]
        Cp_l = data(3);     % Heat Capacity at Constant Pressure
        Cv_l = data(4);     % Heat Capacity at Constant Volume

        % Derivative of specific volume with temperature for constant pressure
        dnu_dT_P = (-1 / (rho_l1^2)) * PropsSI('d(D)/d(T)|P', 'T', T_t, 'P', P1*10^6, 'N2O');

        % Derivative of Compressibility Factor with temperature at constant pressure
        dZ_dT_P = ((P1*10^6)/(R)) * ...
            (T_t^-1 * dnu_dT_P - (rho_l1^-1)/(T_t^2));

        % Derivative of Compressibility Factor with temperature for constant density
        dZ_dT_rho = (1 / (rho_l1 * R)) * (T_t^-1 * ...
            PropsSI('d(P)/d(T)|D', 'T', T_t, 'P', P1*10^6, 'N2O') - (P1*10^6)/T_t^2);

        Z = (P1 * 10^6 * rho_l1^-1) / (R * T_t);   % Compressibility Factor
        gamma = Cp_l / Cv_l;                    % Specific Heats Ratio

        % Isentropic Power Law Exponent for Real Gases
        n = gamma * ((Z + T_t * dZ_dT_rho)/(Z + T_t * dZ_dT_P));

        % Compressibility Correction Factor
        upsilon_line = sqrt(((P1*10^6) / (2 * (P1*10^6 - P2*10^6))) ...
            * ((2*n)/(n - 1)) * (1 - (P1*10^6 - P2*10^6)/(P1*10^6))^(2/n) ...
            * (1 - (1 - (P1*10^6 - P2*10^6)/(P1*10^6))^((n-1)/n)));

        % Mass Flow for Single-Phase Compressible Model [kg/s]
        mdot_SPC = Cd_inj * upsilon_line * A_inj * ...
            sqrt(2 * rho_l1 * (P1*10^6 - P2*10^6));

        % Critical Mass Flow for SPC model [kg/s]
        mdot_SPC_crit = Cd_inj * A_inj * ...
            sqrt(n * rho_l1 * (P1*10^6) * (2 / (n+1))^((n + 1)/(n - 1)));

        crit_ratio = ((2/(n+1))^(n/(n-1)));
        pres_ratio = P2 / P1;
        
        % If pressure ratio exceeds critical value, output critical value
        % instead for mass flow
        if pres_ratio < crit_ratio
            out = mdot_SPC_crit;
        else
            out = mdot_SPC;
        end
    elseif strcmp(model, 'fml')
        
        P1 = P_tank - delta_p;
        P2 = Pc;
        
%         data = N2OSatInterpolator (N2OSat, P2, 'P', {'rho_v','rho_l','s_v','s_l'});
%         
%         rho_2v = data(1);
%         rho_2l = data(2);
%         s_2v = data(3)/MW;
%         s_2l = data(4)/MW;

        data_v = num2cell(PropsSI({'S', 'D'}, 'P', P2, 'Q', 1, 'N2O'));
        data_l = num2cell(PropsSI({'S', 'D'}, 'P', P2, 'Q', 0, 'N2O'));
        
        [s_2v, rho_2v] = data_v{:};
        [s_2l, rho_2l] = data_l{:};
        
        x2 = (s_tank - s_2l)/(s_2v - s_2l);
        
        S = (rho_2l/rho_2v)^(1/3);
        
        alpha_2 = 1 / (1 + ((1-x2)/(x2)) * S * ((rho_2v)/(rho_2l)));
        
        spc = MassFlow (TankData, P_tank, const, s_tank, rho_tank, 'spc');
        aux = MassFlow (TankData, P_tank, const, s_tank, rho_tank, 'hemc');
        hem = aux.mdotc;
        
        out = (1-alpha_2) * spc + alpha_2 * hem;
    else
        % Error message for non-valid model string input
        error('Input string for mass flow model is not valid!');
    end
end

function [x] = secant(fun, x1, in)
    %SECANT is a zero-finding function based on the secant method
    %   'fun' is the function handle for which the zero is desired
    %   'x1' is the initial guess
    %   'in' is a struct containing any additional inputs 'fun' might require
    %   'x' is the value of the zero
    x_eps = x1*0.005; % Set the tolerance to be 0.5% of initial guess
    x2 = x1-x1*0.01; % Set a second point 1% away from the original guess
    F1 = fun(x1, in); % Evaluate function at x1
    F2 = fun(x2, in); % Evaluate function at x2
    kk = 1; % Set up counter
    kk_max = 1000;
    while abs(x2-x1)>=x_eps && kk<kk_max % While error is too large and counter is less than max
        x3 = x2 - F2*(x2-x1)/(F2-F1);
%           if x3 <= 0 %to prevent the secant method from trying to do calculations with negative temperature or Mach
%               x3 = 0.1;
%           end
        x1 = x2; % Move everything forward
        x2 = x3;
        F1 = F2;
        F2 = fun(x2, in);
        kk = kk+1;
    end
    x = x2;
end


function A = Aerror(M, in) % Finds the difference between the estimated and actual nozzle ratio
    A = (1/M^2)*(2./(in.k+1).*(1+(in.k-1)./2.*M.^2)).^((in.k+1)./(in.k-1)) - in.A^2;
end


