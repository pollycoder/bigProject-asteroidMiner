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
coeAsteroid0 = [4.374587943314110e+08, 0.134098123850821, 0.0540505062211469, 2.61854482481308, 4.00216803342331, 3.31673390721605];
coeEarth0 = [1.495484423703440e+08, 0.0163866660358080, 5.40080930104537e-05, 3.71828887427766, 4.38789863130065, 6.20499744208261];
coeMars0 = [2.279254603773820e+08, 0.0934491898618057, 0.0322523881233316, 0.863747331544666, 5.00261081874214, 1.94894057775148];
muMars = 4.282837521400000e+04;
muSun = 1.327124400180000e+11;
day = 86400;
g0 = 9.806650000000000e-3;
Isp = 3000;
RMars = 3.389920000000000e+03;
rpMin = 300 + RMars;

% Unit transform
% To make the calculation faster and preciser.
% From now on, all the calculation will be completed 
% in new unit system.
lUnit = 1 / coeEarth0(1);                                   % Length (AU)
tUnit = 1 / (2 * pi * sqrt(coeEarth0(1) ^ 3 / muSun));      % Time (y)
vUnit = lUnit / tUnit;                                      % Velocity (AU/y)
aUnit = lUnit / tUnit^2;                                    % Acceleration (AU/y^2)
muUnit = lUnit ^ 3 / tUnit ^ 2;                             % Mu (AU^3/y^2)
coeUnit = [lUnit, ones(1, 5)];                              % Change the unit of orbit elements quickly

% New unit - They will be set as global constant
muSunNew = muSun * muUnit;
muMarsNew = muMars * muUnit;
g0New = g0 * aUnit;
IspNew = Isp * tUnit;

coeEarth0New = coeEarth0 .* coeUnit;
coeMars0New = coeMars0 .* coeUnit;
coeAsteroid0New = coeAsteroid0 .* coeUnit;

rpMinNew = rpMin * lUnit;

tWaitUpper = 1825 * day;
tTotalUpper = 5475 * day;
tWaitUpperNew = tWaitUpper * tUnit;
lb = [0, 3, 5, 7, 7, 10, 10, 0, 0, 0, 0, 0, 0]';
ub = [3, 15, 15, 15, 15, 15, 15, 1, 1, 1, 2 * pi, 2 * pi, 2 * pi]';

%% Optimize fuel - global - PSO
%options = optimoptions("particleswarm", "SwarmSize", 1000, ...
%                       'UseParallel', true, 'MaxIterations', 1000, ...
%                       'HybridFcn', 'patternsearch', 'Display', 'iter');
%[X, init_result, exitflag] = particleswarm(@GA4_obj, 13, lb, ub, options);


%% Optimize fuel - global -GA
options = optimoptions("ga", "ConstraintTolerance", 1e-10, "CreationFcn", ...
                       "gacreationlinearfeasible", "CrossoverFcn", "crossoverlaplace", ...
                       "Display", "iter", "HybridFcn", "patternsearch", 'UseParallel', true);
[X,fval,exitflag,output,population,~] = ga(@GA4_obj, 13, [], [], [], [], lb, ub, [], [], options);


%% Optimize fuel - local
options = optimset('MaxIter', 10000);
[X, result] = fminsearch(@GA4_obj, X, options);
fprintf("J=%f\n",result);
