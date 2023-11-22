clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Big Project - Asteroid Mining
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constant
% Unit:
% Length: km
% Angle: rad
% Time: s
load("data/constant.mat");
load("data/coeAsteroid0.mat");
load("data/coeEarth0.mat");
load("data/coeMars0.mat");

%% Unit transform
% To make the calculation faster and preciser.
% From now on, all the calculation will be completed 
% in new unit system.
lUnit = coeEarth0(1);                               % Length (AU)
tUnit = 2 * pi * sqrt(coeEarth0(1) ^ 3 / muSun);    % Time (y)
vUnit = lUnit / tUnit;                              % Velocity (AU/d)
muUnit = lUnit ^ 3 / tUnit ^ 2;                     % Mu (AU^3/d^2)
coeUnit = [lUnit, ones(1, 5)];                      % Change the unit of orbit elements quickly

% New unit - They will be set as global constant
global muSunNew muMarsNew coeEarth0New coeMars0New coeAsteroid0New
muSunNew = muSun / muUnit;
muMarsNew = muMars / muUnit;

coeEarth0New = coeEarth0 ./ coeUnit;
coeMars0New = coeMars0 ./ coeUnit;
coeAsteroid0New = coeAsteroid0 ./ coeUnit;

tWaitUpperNew = tWaitUpper / tUnit;
tTotalUpperNew = tTotalUpper / tUnit;

