%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first solution: mono-GA method
% One GA for departure, one GA for return.
% Caution: index:
% - New: new unit
% - Int: international unit
% - Km: especially for length, velocity and mu, 
%       because the length unit here is km
% X(1)~X(5): Time
% X(6)~X(7): rp
% X(8)~X(9): Phi
% X(10): mf
% X: new unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = monoGA_obj(X)
tol = 1e-12;

% Penalty for rp
hpmin = evalin("base", 'rpmin');
RMars = evalin("base", 'RMars');
rpMin = hpmin + RMars;
if X(6) < rpMin || X(7) < rpMin
    J = 1;
    return
end

% Transform the variables from main workspace to here
coeAsteroid0 = evalin("base", 'coeAsteroid0');
coeEarth0 = evalin("base", 'coeEarth0');
coeMars0 = evalin("base", 'coeMars0');
muSun = evalin("base", 'muSun');
muMars = evalin("base", 'muMars');
Isp = evalin("base", 'Isp');

% Unit transform
% To make the calculation faster and preciser.
% From now on, all the calculation will be completed 
% in new unit system.
lUnit = 1 / coeEarth0(1);                                   % Length (AU)
tUnit = 1/ (2 * pi * sqrt(coeEarth0(1) ^ 3 / muSun));       % Time (y)
vUnit = lUnit / tUnit;                                      % Velocity (AU/d)
muUnit = lUnit ^ 3 / tUnit ^ 2;                             % Mu (AU^3/d^2)
coeUnit = [lUnit, ones(1, 5)];                              % Change the unit of orbit elements quickly

% New unit - They will be set as global constant
muSunNew = muSun * muUnit;
muMarsNew = muMars * muUnit;

coeEarth0New = coeEarth0 .* coeUnit;
coeMars0New = coeMars0 .* coeUnit;
coeAsteroid0New = coeAsteroid0 .* coeUnit;

tWaitUpperNew = tWaitUpper * tUnit;
tTotalUpperNew = tTotalUpper * tUnit;

% Initial mass
mDry = 500;                                                 % Initial dry mass (kg)                     
mFuel = 500;                                                % Initial fuel mass (kg)
mTotal0 = mDry + mFuel;                                     % Initial total mass (kg)

% Departure: Earth -> Mars
tEMNew = X(2) - X(1);                                       % Transfer time
[rE0New, vE0New] = coe2rv(coeEarth0New, muSunNew, tol);     % RV of Earth (t=0)
[rM0New, vM0New] = coe2rv(coeMars0New, muSunNew, tol);      % RV of Mars (t=0)

[rEt0New, vEt0New] = rv02rvf(rE0New, vE0New, ...
                             X(1), muSunNew);               % RV of Earth (t=t0)
[rMt1New, vMt1New] = rv02rvf(rM0New, vM0New, ...
                             X(2), muSunNew);               % RV of Mars (t=t1)

[vt0New, vt11New] = LambSol(rEt0New, rMt1New, ...
                            tEMNew, muSunNew);              % Lambert problem 1: E->M

dv1New = vt0New - vEt0New;                                  % First impulse (t=t0)
dv1 = norm(dv1New) / vUnit;                                 % Unit transform to fit the impulse solver
[mTotalt0, dmt0] = impulseFuel(mTotal0, dv1, Isp);          % Mass change (t=t0)
mFuel = mFuel - dmt0;                                       % Fuel cost (t=t0)

if mFuel < 0                                                % Penalty, to stop calculation in time if condition is not satisfied
    J = 1;                                                  % J < 0, therefore penalty > 0
    return
end

% GA-1: SOI
% vt11 would become the velocity for GA
[vt12New, ~] = SOI_after(vt11New, vMt1New, muMarsNew, ...
                         X(6), X(8));                       % SOI

% Arrival: Mars->Asteroid
tMA = X(2) - X(1);
[rA0New, vA0New] = coe2rv(coeAsteroid0New, muSunNew, tol);  % RV of Asteroid (t=0)
[rAt2New, vAt2New] = rv02rvf(rA0New, vA0New, ...
                             X(2), muSunNew);               % RV of Asteroid (t=t2)

[vt13New, vt2New] = LambSol(rMt1New, rAt2New, ...
                            tMA, muSunNew);                 % Lambert problem 2: M->A

dv2New = vt13New - vt12New;                                 % Second impulse (t=t1)
dv2 = norm(dv2New) / vUnit;                                 % Unit transform to fit the impulse solver
[mTotalt1, dmt1] = impulseFuel(mTotalt0, dv2, Isp);         % Mass change (t=t1)
mFuel = mFuel - dmt1;                                       % Fuel cost (t=t1)

dv3New = vAt2New - vt2New;                                  % Third impulse (t=t2)
dv3 = norm(dv3New) / vUnit;                                 % Unit transform to fit the impulse solver
[mTotalt21, dmt2] = impulseFuel(mTotalt1, dv3, Isp);        % Mass change (t=t2)
mFuel = mFuel - dmt2;                                       % Fuel cost (t=t2)

if mFuel < 0                                                % Penalty, to stop calculation in time if condition is not satisfied
    J = 1;                                                  % J < 0, therefore penalty > 0
    return
end

% Sampling
mTotalt22 = mTotalt21 + X(10);                              % Add sample mass
mDry = mDry + X(10);                                        % Sample mass is included in dry mass

% Return: A->M



J = -X(10);
end