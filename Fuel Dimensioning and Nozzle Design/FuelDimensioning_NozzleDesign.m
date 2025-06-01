%% Fuel Dimensions and Nozzle Design Code

%% Intro

% Thesis code developed by Guilherme Folgado, 92685

% This script allows dimensioning of a BATES style fuel grain from several
% preliminary performance parameters that can be computed using other
% scripts developed by the author (ie HybridPrelimAnalysisv2); also
% computes critical nozzle pressure ratios to be used in the nozzle flow
% model (preliminary version of said model is implemented too).

clc
clear all

%% Initialize variables

data = [];
files = {};

% Set both fuels and oxidizer. Variables are name, mass fraction, temp,
% enthalpy of formation and exploded chemical formula.
fuels = {'C32H66(a)', 0.9, 298.15, 0, ''; 'stacid', 0.1, 298.15, -912, 'C 18 H 36 O 2'};
oxs = {'N2O(L),298.15K', 1, 298.15, 0, ''};

OF = 8;                 % Oxidizer to Fuel ratio
Pc = 30;                % CC pressure [bar] 
Patm = 0.9197;          % Atmospheric pressure for optimized nozzle [bar]
F = 10000;               % Expected Thrust [N]
comb_eff = .95;         % Combustion efficiency
rho_f = 920;            % Fuel Density [kg/m^3]
a = 0.155 * 10^-3;      % Regression rate coefficient
n = 0.5;              % Regression rate exponent
t_b = 12;                % Burn time [s]
wa = .005;              % Aditional web thickness [m]
nozzle_ratio = 5;   % Nozzle throat to exit area ratio

% Generate data tables from user inputs above, to then use for computing
% fuel grain dimensions, and critical pressure ratios

[data, files] = CEA_TableGenerator(Pc,OF,nozzle_ratio,fuels,oxs);

% Get these values from CEA prelim analysis code

c_star = data(15);           % Characteristic velocity [m/s]
Cf = data(17);              % Thrust Coefficient

% REMINDER: these output values from CEA still need to have correction
% factors applied to them; see pages 43-44 and 60-61 of document
% Development of a Hybrid Sounding Rocket Motor by Genevieve

%% Fuel Dimensioning Results

A_t = F / (Cf * Pc * 10^5);                             % Nozzle throat area [m^2]
D_t = 2 * sqrt(A_t / pi);                               % Nozzle throat diameter [m]
r_t = D_t / 2;
A_e = A_t * nozzle_ratio;                               % Nozzle exit area [m^2]
D_e = 2 * sqrt(A_e / pi);                               % Nozzle exit diameter [m]
r_e = D_e / 2;
mdot_nozzle = (Pc * 10^5 * A_t) / (comb_eff * c_star);  % Nozzle mass flow [kg/s]
mdot_fuel = mdot_nozzle / (OF + 1);                     % Fuel mass flow [kg/s]
mdot_ox = OF * mdot_fuel;                               % Oxidizer mass flow [kg/s]

I_t = F * t_b;                                          % Total Impulse [Ns]
m_fuel_t = mdot_fuel * t_b;                             % Total initial fuel mass (without additional web thickness) [kg]
m_ox_t = mdot_ox * t_b;                                 % Total initial ox mass (rough estimate)
V_fuel_i = m_fuel_t / rho_f;                            % Respective fuel volume [m^3]

% Function for regression rate to be integrated

dr_f_dt = @(t, r_f) a * (mdot_ox ./ (pi * r_f .^2)) .^ n;

r_f0 = 0.03;        % Initial value for fuel port radius [m]

tspan = [0, t_b];   % Time interval for integration, from 0 until burn time

% Solve differential equation

[t, r_f] = ode45(dr_f_dt, tspan, r_f0);

% Integrate differential eq solution to get theoretical web thickness

wb = r_f(end) - r_f(1);

wt = wb + wa;                           % Total web thickness [m]
Rf = wt + r_f0;                         % Initial grain outer radius [m]
Df = 2 * Rf;
d_f0 = 2 * r_f0;

Lg = (4 * V_fuel_i) / (pi() * (Df^2 - d_f0^2));   % Fuel grain length [m]

port_throat_ratio = (pi() * r_f0^2) / (A_t);