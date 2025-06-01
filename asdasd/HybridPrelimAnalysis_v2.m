%% Hybrid Rocket Engine Prelim Analysis Code

%% Intro

% Thesis Work Done by Guilherme Folgado, IST student number 92685 

% This is the preliminary design tool, to allow the user to 
% evaluate important operating parameters (such as nozzle expansion ratio, 
% characteristic velocity, and specific impulse) from user inputs (chamber 
% pressure, optimal expansion altitude, OF ratios, fuels and oxidizers used).
% The computations are then performed, and the resulting data is presented 
% in graphs.
% Also allows just viewing data from provided .out files, so long as naming
% convention is followed ('pressxx.out', where xx is the Pc value) and the
% analysis was performed for equilibrium conditions (infinite area CC)

%% Select New Analysis or Run Previous Analysis

clc;
clear all;
option = input('Do you wish to run a new analysis (N) or to load the previous analysis (L)? ', 's');

switch option
%% User Inputs
% Gets inputs from user for analysis (chamber pressures, OF ratios,
% altitude to optimize nozzle expansion, fuels and oxidizers to be used);
% returns those inputs and also the list of reagent species names included
% in the CEA thermo library

    case 'N'
        
        [Pc, OF, Patm, fuels, oxs, reagents] = UserInputs();

%% Running CEA, saving data
% From inputs, creates input files for CEA application, automates running
% each .inp file in CEA, and data collection for performance parameters,
% that are stored in the prelim_data matrix; FileList is a list with all
% the created .inp and .out files from CEA

        [FileList, prelim_data] = CEAPrelimMaster(Pc, OF, Patm, fuels, oxs);

%% Plotting
% From CEA analysis results, queries user as to which performance
% parameters should be plotted in graphs for visualization

        Plotting(prelim_data, Pc);

%% Load previous Analysis
% Load all the .out files in the CEA folder with naming convention of
% 'pressxx.out' where xx is the chamber pressure value. Allows the user to
% bypass inputting variables, and also to just provide .out files for
% displaying analysis results. NOTE: only use .out files stemming from
% analysis using equilibrium conditions, or data may not load properly

    case 'L'
        
        [FileList, Pc] = getfiles();
        prelim_data = PrelimData(FileList);
        Plotting(prelim_data, Pc)
        
%% Invalid input (neither new or load selected)        

    otherwise
        
        fprintf('Invalid input. Must be either N or L.\n');
        
end