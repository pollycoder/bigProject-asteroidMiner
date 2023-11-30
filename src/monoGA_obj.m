%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first solution: mono-GA method
% One GA for departure, one GA for return.
% Caution: index:
% - New: new unit
% - Int: international unit
% - Km: especially for length, velocity and mu, 
%       because the length unit here is km
% X(1)~X(6): Time
% X(7)~X(8): rp
% X(9)~X(10): Phi
% X(11): mf
% X: new unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = monoGA_obj(X)
tol = 1e-12;
penalty = 1e20;                        % Since the result should be negative, any positive number could be penalty
% Penalty for time 
if ~all(X(:) > 0)
    J = penalty;
    return
end
dX = diff(X);
if ~all(dX(1:5) > 0)
    %warning("时间出现反转。Penalty.");
    J = penalty;
    return;
end

if X(1) > 4.9990 || X(6) > 14.9990
    %warning("任务时间超过，Penalty.");
    J = penalty;
    return
end

%% Constant
% Unit:
% Length: km
% Angle: rad
% Time: s
coeAsteroid0 = [4.374587943314110e+08, ...
                0.134098123850821, ...
                0.0540505062211469, ...
                2.61854482481308, ...
                4.00216803342331, ...
                3.31673390721605];
coeEarth0 = [1.495484423703440e+08, ...
             0.0163866660358080, ...
             5.40080930104537e-05, ...
             3.71828887427766, ...
             4.38789863130065, ...
             6.20499744208261];
coeMars0 = [2.279254603773820e+08, ...
            0.0934491898618057, ...
            0.0322523881233316, ...
            0.863747331544666, ...
            5.00261081874214, ...
            1.94894057775148];
muMars = 4.282837521400000e+04;
muSun = 1.327124400180000e+11;
Isp = 3000;
g0 = 9.806650000000000e-3;
RMars = 3.389920000000000e+03;
hpmin = 300;


% Unit transform
% To make the calculation faster and preciser.
% From now on, all the calculation will be completed 
% in new unit system.
lUnit = 1 / coeEarth0(1);                                   % Length (AU)
tUnit = 1 / (2 * pi * sqrt(coeEarth0(1) ^ 3 / muSun));      % Time (y)
vUnit = lUnit / tUnit;                                      % Velocity (AU/y)
aUnit = lUnit / tUnit^2;
muUnit = lUnit ^ 3 / tUnit ^ 2;                             % Mu (AU^3/y^2)
coeUnit = [lUnit, ones(1, 5)];                              % Change the unit of orbit elements quickly

% New unit - They will be set as global constant
muSunNew = muSun * muUnit;
muMarsNew = muMars * muUnit;
IspNew = Isp * tUnit;
g0New = g0 * aUnit;

coeEarth0New = coeEarth0 .* coeUnit;
coeMars0New = coeMars0 .* coeUnit;
coeAsteroid0New = coeAsteroid0 .* coeUnit;

% Penalty for rp - new unit
rpMin = hpmin + RMars;
rpMinNew = rpMin * lUnit;
if X(7) < rpMinNew || X(8) < rpMinNew
    %waring("借力高度不够。Penalty.")
    J = penalty;
    return
end

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
    dvt0Norm = 0;
    dvt0NormNew = 0;
else
    dvt0Norm = dvt0Norm - 4;
    dvt0NormNew = dvt0Norm * vUnit;
end
[mTotalt0, dmt0] = impulseFuel(mTotal0, dvt0NormNew, ...
                               IspNew, g0New);              % Mass change (t=t0)
mFuel = mFuel - dmt0;                                       % Fuel cost (t=t0)


% GA-1: SOI (t1)
% vt11 would become the velocity for GA
[vt12New, ~] = SOI_after(vt11New, vMt1New, muMarsNew, ...
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


% GA-1: SOI (t1)
% vt11 would become the velocity for GA
[vt42New, ~] = SOI_after(vt41New, vMt4New, muMarsNew, ...
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
    dvt5Norm = 0;
    dvt5NormNew = 0;
else
    dvt5Norm = dvt0Norm - 4;
    dvt5NormNew = dvt0Norm * vUnit;
end
[mTotalt5, dmt5] = impulseFuel(mTotalt4, dvt5NormNew, ...
                               IspNew, g0New);              % Mass change (t=t5)
mFuel = mFuel - dmt5;

%{
if mFuel < 0                                                % Penalty, to stop calculation in time if condition is not satisfied
    %warning("脉冲6,燃料耗尽。Penalty.");
    J = penalty;                                            % J < 0, therefore penalty > 0
    return
end
%}

J = -mFuel;

end