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
lb = [0, 2, 5, 6, 9, 9, 0, 0, 0, 0, 0]';
ub = [2, 5, 9, 15, 15, 15, 10, 10, 2 * pi, 2 * pi, 10]';

%%
options = optimoptions("particleswarm", "SwarmSize", 1000, 'UseParallel', true, 'MaxIterations', 1000, 'Display', 'iter');
[X, init_result, exitflag] = particleswarm(@monoGA_obj1, 11, lb, ub, options);

%%
options = optimset('MaxIter', 10000, 'Display', 'iter');
[X, result] = fminsearch(@monoGA_obj1, X, options);

%%
X(11) = 895.79;
options = optimset('MaxIter', 10000, 'Display', 'iter');
[X, result] = fminsearch(@monoGA_obj, X, options);
fprintf("m=%f\n", X(11));
fprintf("J=%f\n",result);


%%
tol = 1e-20;
% Initial mass
mDry = 500;                                                 % Initial dry mass (kg)                     
mFuel = 500;                                                % Initial fuel mass (kg)
mTotal0 = mDry + mFuel;                                     % Initial total mass (kg)

% Departure: E -> M (t0-t1)
tEMNew = X(2) - X(1);                                       % Transfer time
[rE0New, vE0New] = coe2rv(coeEarth0New, muSunNew, tol);     % RV of Earth (t=0)
[rM0New, vM0New] = coe2rv(coeMars0New, muSunNew, tol);      % RV of Mars (t=0)

[rEt0New, vEt0New] = rv02rvf(rE0New, vE0New, ...
                             X(1), muSunNew);               % RV of Earth (t=t0)
[rMt1New, vMt1New] = rv02rvf(rM0New, vM0New, ...
                             X(2), muSunNew);               % RV of Mars (t=t1)

[vt0New, vt11New] = LambSol(rEt0New, rMt1New, ...
                            tEMNew, muSunNew);              % Lambert problem 1: E->M

dvt0New = vt0New - vEt0New;                                 % 1st impulse (t=t0)
dvt0NormNew = norm(dvt0New);                                % Unit transform to fit the impulse solver
dvt0Norm = dvt0NormNew / vUnit;
if dvt0Norm < 4
    dvt0New = zeros(3, 1);
    dvt0Norm = 0;
    dvt0NormNew = 0;
else
    dvt0Vector = dvt0New / dvt0NormNew;
    dvt0New = dvt0New - 4 * dvt0Vector * vUnit;
    dvt0Norm = dvt0Norm - 4;
    dvt0NormNew = dvt0Norm * vUnit;
end
[mTotalt0, dmt0] = impulseFuel(mTotal0, dvt0NormNew, ...
                               IspNew, g0New);              % Mass change (t=t0)
mFuel = mFuel - dmt0;                                       % Fuel cost (t=t0)


% GA-1: SOI (t1)
% vt11 would become the velocity for GA
[vt12New, dvt1GANew] = SOI_after(vt11New, vMt1New, muMarsNew, ...
                         X(7), X(9));                       % SOI

% Arrival: M->A (t1-t2)
tMA = X(3) - X(2);
[rA0New, vA0New] = coe2rv(coeAsteroid0New, muSunNew, tol);  % RV of Asteroid (t=0)
[rAt2New, vAt2New] = rv02rvf(rA0New, vA0New, ...
                             X(3), muSunNew);               % RV of Asteroid (t=t2)

[vt13New, vt2New] = LambSol(rMt1New, rAt2New, ...
                            tMA, muSunNew);                 % Lambert problem 2: M->A

dvt1New = vt13New - vt12New;                                % 2nd impulse (t=t1)
dvt1NormNew = norm(dvt1New);                                % Unit transform to fit the impulse solver
[mTotalt1, dmt1] = impulseFuel(mTotalt0, dvt1NormNew, ...
                               IspNew, g0New);              % Mass change (t=t1)
mFuel = mFuel - dmt1;                                       % Fuel cost (t=t1)

dvt2New = vAt2New - vt2New;                                 % 3rd impulse (t=t2)
dvt2NormNew = norm(dvt2New);                                % Unit transform to fit the impulse solver
[mTotalt21, dmt2] = impulseFuel(mTotalt1, dvt2NormNew, ...
                                IspNew, g0New);             % Mass change (t=t2)
mFuel = mFuel - dmt2;                                       % Fuel cost (t=t2)


% Sampling (t2)
mTotalt22 = mTotalt21 + X(11);                              % Add sample mass
mDry = mDry + X(11);                                        % Sample mass is included in dry mass

% Return: A->M (t3-t4)
tAM = X(5) - X(4);
[rAt3New, vAt3New] = rv02rvf(rA0New, vA0New, ...
                             X(4), muSunNew);               % RV of Asteroid (t=t3)
[rMt4New, vMt4New] = rv02rvf(rM0New, vM0New, ...
                             X(5), muSunNew);               % RV of Mars (t=t4)

[vt3New, vt41New] = LambSol(rAt3New, rMt4New, ...
                            tAM, muSunNew);                 % Lambert problem 3: A->M

dvt3New = vt3New - vAt3New;                                 % 4th impulse (t=t3)
dvt3NormNew = norm(dvt3New);                                % Unit transform to fit the impulse solver
[mTotalt3, dmt3] = impulseFuel(mTotalt22, dvt3NormNew, ...
                               IspNew, g0New);              % Mass change (t=t3)
mFuel = mFuel - dmt3;                                       % Fuel cost (t=t3)


% GA-1: SOI (t4)
% vt11 would become the velocity for GA
[vt42New, dvt4GANew] = SOI_after(vt41New, vMt4New, muMarsNew, ...
                         X(8), X(10));                      % SOI

% Return: M->E (t4-t5)
tME = X(6) - X(5);
[rEt5New, vEt5New] = rv02rvf(rE0New, vE0New, ...
                             X(6), muSunNew);               % RV of Earth (t=t5)

[vt43New, vt5New] = LambSol(rMt4New, rEt5New, ...
                            tME, muSunNew);                 % Lambert problem 4: M->E

dvt4New = vt43New - vt42New;                                % 5th impulse (t=t4)
dvt4NormNew = norm(dvt4New);                                % Unit transform to fit the impulse solver
[mTotalt4, dmt4] = impulseFuel(mTotalt3, dvt4NormNew, ...
                               IspNew, g0New);              % Mass change (t=t4)
mFuel = mFuel - dmt4;

dvt5New = vEt5New - vt5New;                                 % 6th impulse (t=t5)
dvt5NormNew = norm(dvt5New);                                % Unit transform to fit the impulse solver
dvt5Norm = dvt5NormNew / vUnit;
if dvt5Norm < 4
    dvt5New = zeros(3, 1);
    dvt5Norm = 0;
    dvt5NormNew = 0;
else
    dvt5Vector = dvt5New / dvt5NormNew;
    dvt5New = dvt5New - 4 * dvt5Vector * vUnit;
    dvt5Norm = dvt5Norm - 4;
    dvt5NormNew = dvt5Norm * vUnit;
end
[mTotalt5, dmt5] = impulseFuel(mTotalt4, dvt5NormNew, ...
                               IspNew, g0New);              % Mass change (t=t5)
mFuel = mFuel - dmt5;

dvt0 = dvt0New / vUnit;
dvt1 = dvt1New / vUnit;
dvt2 = dvt2New / vUnit;
dvt3 = dvt3New / vUnit;
dvt4 = dvt4New / vUnit;
dvt5 = dvt5New / vUnit;
dvt1GA = dvt1GANew / vUnit;
dvt4GA = dvt4GANew / vUnit;

vA0 = vA0New / vUnit;
vAt2 = vAt2New / vUnit;
vAt3 = vAt3New / vUnit;
vEt0 = vEt0New / vUnit;
vEt5 = vEt5New / vUnit;
vMt1 = vMt1New / vUnit;
vMt4 = vMt4New / vUnit;
vt0 = vt0New / vUnit;
vt11 = vt11New / vUnit;
vt12 = vt12New / vUnit;
vt13 = vt13New / vUnit;
vt2 = vt2New / vUnit;
vt3 = vt3New / vUnit;
vt41 = vt41New / vUnit;
vt42 = vt42New / vUnit;
vt43 = vt43New / vUnit;
vt5 = vt5New / vUnit;

rA0 = rA0New / lUnit;
rAt2 = rAt2New / lUnit;
rAt3 = rAt3New / lUnit;
rE0 = rE0New / lUnit;
rEt0 = rEt0New / lUnit;
rEt5 = rEt5New / lUnit;
rM0 = rM0New / lUnit;
rMt1 = rMt1New / lUnit;
rMt4 = rMt4New / lUnit;

X_int = X;
X_int(1:6) = X(1:6) / tUnit / day;
X_int(9) = mod(X_int(9), 2 * pi);
X_int(10) = mod(X_int(10), 2 * pi);


%% Plot
% Earth, Mars, Asteroid
n=100;
f = linspace(0, 2 * pi, n);
for i = 1:length(f)
    coeENew = coeEarth0New;
    coeENew(6) = f(i);
    rE(:, i) = coe2rv(coeENew, muSunNew, tol);
end
for i = 1:length(f)
    coeMNew = coeMars0New;
    coeMNew(6) = f(i);
    rM(:, i) = coe2rv(coeMNew, muSunNew, tol);
end
for i = 1:length(f)
    coeANew = coeAsteroid0New;
    coeANew(6) = f(i);
    rA(:, i) = coe2rv(coeANew, muSunNew, tol);
end
plot3(rE(1,:), rE(2,:), rE(3,:), 'LineWidth', 1.5, 'Color', 'k');hold on
plot3(rM(1,:), rM(2,:), rM(3,:), 'LineWidth', 1.5, 'Color', 'k');hold on
plot3(rA(1,:), rA(2,:), rA(3,:), 'LineWidth', 1.5, 'Color', 'k');hold on


% Trajectory
r0 = rEt0New;
v0 = vt0New;
n1=1000;
t01 = linspace(0, X(2) - X(1), n1);
for i=1:length(t01)
    [r(:, i), ~] = rv02rvf(r0, v0, t01(i), muSunNew);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'r', 'LineStyle','--');hold on
plot3(r(1,1), r(2,1), r(3,1),'g*','LineWidth', 2);hold on
plot3(r(1, end), r(2,end), r(3,end),'b*','LineWidth', 2);hold on
text(r(1,1), r(2,1), r(3,1), 'Departure');
text(r(1, end), r(2,end), r(3,end), 'GA-Mars-1');

r0 = rMt1New;
v0 = vt13New;
t12 = linspace(0, X(3) - X(2), n1);
for i=1:length(t12)
    [r(:, i), ~] = rv02rvf(r0, v0, t12(i), muSunNew);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'm', 'LineStyle','--');hold on
plot3(r(1, end), r(2,end), r(3,end),'m*','LineWidth', 2);hold on
text(r(1, end), r(2,end), r(3,end), 'Arrival-Pysche-Mining');

r0 = rAt2New;
v0 = vAt2New;
t23 = linspace(0, X(4)- X(3), n1);
for i=1:length(t12)
    [r(:, i), ~] = rv02rvf(r0, v0, t23(i), muSunNew);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'g', 'LineStyle','--');hold on
plot3(r(1, end), r(2,end), r(3,end),'g*','LineWidth', 2);hold on
text(r(1, end), r(2,end), r(3,end), 'Return-Psyche');

r0 = rAt3New;
v0 = vt3New;
t34 = linspace(0, X(5)- X(4), n1);
for i=1:length(t12)
    [r(:, i), ~] = rv02rvf(r0, v0, t34(i), muSunNew);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'c', 'LineStyle','--');hold on
plot3(r(1, end), r(2,end), r(3,end),'b*','LineWidth', 2);hold on
text(r(1, end), r(2,end), r(3,end),'GA-Mars-2');

r0 = rMt4New;
v0 = vt43New;
t45 = linspace(0, X(6)- X(5), n1);
for i=1:length(t12)
    [r(:, i), ~] = rv02rvf(r0, v0, t45(i), muSunNew);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'b', 'LineStyle','--');hold on
plot3(r(1, end), r(2,end), r(3,end),'r*','LineWidth', 2);hold on 
text(r(1, end), r(2,end), r(3,end),'Arrival-Earth');

plot3(0,0,0,'k*','LineWidth', 3);
text(0,0,0,'Sun');

axis equal