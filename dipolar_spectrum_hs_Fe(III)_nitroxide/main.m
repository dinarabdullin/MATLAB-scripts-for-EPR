% Simulation of a dipolar spectrum for a high-spin iron(III)/nitroxide pair
% High-spin iron(III) is treated as a spin-5/2.  The problem is solved via 
% diagonalization of the total Hamiltonian.

% Reset MATLAB memory and command line
clear all
clc

%% Constants
Const.Hz2MHz = 1e-6;
Const.MHz2Hz = 1e6;
Const.GHz2MHz = 1e3;
Const.mT2T = 1e-3;
Const.T2mT = 1e3;
Const.nm2m = 1e-9;
Const.deg2rad = pi / 180;
Const.rad2deg = 180 / pi;
Const.wn2MHz = 29979;
Const.plankConstant = 6.626070040e-34; % in J*s
Const.bohrMagneton = 9.274009994e-24; % in J/T
Const.bolzmannConstant = 1.38064852e-23; % in J/K
Const.ge = 2.0023; % free electron g-factor
Const.vacuumPermeabilityOver4Pi = 1e-7; % in T*m/A
Const.Fez = Const.Hz2MHz * Const.bohrMagneton / Const.plankConstant; % in MHz/T
Const.Fdd = Const.Hz2MHz * Const.vacuumPermeabilityOver4Pi * Const.bohrMagneton^2 / ...
            (Const.plankConstant * Const.nm2m^3); % in MHz

%% Experimental parameters
Exp.mwFreq = 33.700000; % in GHz
Exp.magnField = 1.2025; % in T
Exp.temperature = 300; % in K

%% Parameters of the spin system
% High-spin iron(III)
Spin1.S = 5/2;
Spin1.gActual = [2.0023 2.0023 2.0023];
Spin1.gFrame = [0 0 0] * Const.deg2rad;
Spin1.D = [10 0] * Const.wn2MHz;
Spin1.DFrame = [0 0 0] * Const.deg2rad;
% Calculate effective g-tensor of the high-spin iron(III)
Spin1.g = gTensorIron(Spin1, Exp, Const);
fprintf('Principal g values of the high-spin iron(III):\n');
fprintf('gxx = %.4f, gyy = %.4f, gzz = %.4f\n', Spin1.g(1), Spin1.g(2), Spin1.g(3));
% Nitroxide
Spin2.S = 1/2;
Spin2.g = [2.0023 2.0023 2.0023];
Spin2.gFrame = [0 0 0] * Const.deg2rad;

%% Simulation parameters
Sim.r = 2.50; % in nm
Sim.rStd = 0.0; % in nm
Sim.rFrame = [0 90 0] * Const.deg2rad;
Sim.nSamples = 1e7;
Sim.freq = -40:0.1:40;
freq0 = Const.Fdd * Const.ge * Const.ge / Sim.r^3;
Sim.relFreq = Sim.freq / freq0;

%% Output parameters
Output.freqLim = [-40, 40, 10];
Output.relFreqLim = [-8, 8, 1];
Output.spcLim = [0, 1.05, 0.5];
Output.saveFigures = true;
Output.directory = '';
if (Output.saveFigures)
    % Create a folder to save data
    directory = strcat(pwd, datestr(now,'/yy-mm-dd_HH-MM-SS/'));
	[status, msg, msgID] = mkdir(directory);
    Output.directory = directory;
end

%% Simulation
display('Simulation of the dipolar spectrum...')
% Calculate the dipolar spectrum
spc = calculateSpectrum3(Spin1, Spin2, Exp, Sim, Const);
% Plot the dipolar spectrum vs absolute frequency
plotSpectrum(Sim.freq, spc, Output.freqLim, Output.spcLim, false, ...
             Output.saveFigures, strcat(Output.directory, 'spc.png'));

% Plot the dipolar spectrum vs relative frequency
plotSpectrum(Sim.relFreq, spc, Output.relFreqLim, Output.spcLim, true, ...
             Output.saveFigures, strcat(Output.directory, 'spc2.png'));
display('Done!')


%% Save
M = zeros(size(Sim.freq,2),2); 
M(:,1) = Sim.freq; 
M(:,2) = spc'; 
filename = strcat(Output.directory, 'spc.dat'); 
dlmwrite(filename,M,'delimiter', '\t');