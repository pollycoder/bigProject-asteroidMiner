%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second solution: 4-GA method
% Two GAs for departure, two GAs for return.
% Caution: index:
% - New: new unit
% - Int: international unit
% - Km: especially for length, velocity and mu, 
%       because the length unit here is km
% X(1)~X(8): Time
% X(9)~X(12): rp
% X(13)~X(16): Phi
% X: new unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = GA4_obj(X)
tol = 1e-20;
penalty = 1e20;                       

% Time constraint
if ~all(X(:) > 0)
    J = penalty;
    return
end
dX = diff(X);
if ~all(dX(1:7) > 0)
    J = penalty;
    return
end
if X(1) > 4.9990 || X(6) > 14.9990
    J = penalty;
    return
end

% Constant
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
muEarth = 398600;
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
muEarthNew = muEarth * muUnit;
IspNew = Isp * tUnit;
g0New = g0 * aUnit;

coeEarth0New = coeEarth0 .* coeUnit;
coeMars0New = coeMars0 .* coeUnit;
coeAsteroid0New = coeAsteroid0 .* coeUnit;

% Penalty for rp - new unit
rpMin = hpmin + RMars;
rpMinNew = rpMin * lUnit;
if X(9) < rpMinNew || X(10) < rpMinNew || X(11) < rpMinNew || X(12) < rpMinNew
    J = penalty;
    return
end

% Initial mass
mDry = 500;                                                 % Initial dry mass (kg)                     
mFuel = 500;                                                % Initial fuel mass (kg)
mTotal0 = mDry + mFuel;                                     % Initial total mass (kg)

% t=0
[rE0New, vE0New] = coe2rv(coeEarth0New, muSunNew, tol);     % RV of Earth (t=0)
[rM0New, vM0New] = coe2rv(coeMars0New, muSunNew, tol);      % RV of Mars (t=0)
[rA0New, vA0New] = coe2rv(coeAsteroid0New, muSunNew, tol);  % RV of Asteroid (t=0)

% Departure: E-->E (t0-t1)
tEMNew = X(2) - X(1);
[rEt0New, vEt0New] = rv02rvf(rE0New, vE0New, ...
                             X(1), muSunNew);               % RV of Earth (t=t0)
[rEt1New, vEt1New] = rv02rvf(rE0New, vE0New, ...
                             X(2), muSunNew);               % RV of Mars (t=t1)
[vt0New, vt11New] = LambSol(rEt0New, rEt1New, ...
                            tEMNew, muSunNew);              % Lambert problem 1: E->M

% GA-Transfer: E-->M (t1-t2)
tME1New = X(3) - X(2);
[rMt2New, vMt2New] = rv02rvf(rM0New, vM0New, ...
                             X(3), muSunNew);               % RV of Mars (t=t2)
[vt13New, vt21New] = LambSol(rEt1New, rMt2New, ...
                             tME1New, muSunNew);            % Lambert problem 2: M->M

% Arrival: M-->A (t2-t3)
tEANew = X(4) - X(3);
[rAt3New, vAt3New] = rv02rvf(rA0New, vA0New, ...
                             X(4), muSunNew);               % RV of Psyche (t=t3)
[vt23New, vt3New] = LambSol(rMt2New, rAt3New, ...
                            tEANew, muSunNew);              % Lambert problem 3: M->A

% GA-1: SOI (t1)
[vt121New, vt122New, ~] = SOI_opt(vt11New, vt13New, vEt1New, muEarthNew, X(9), X(13));

% GA-2: SOI (t2)
[vt221New, vt222New, ~] = SOI_opt(vt21New, vt23New, vMt2New, muMarsNew, X(10), X(14));

% Return: A-->M (t4-t5)
tAENew = X(6) - X(5);
[rAt4New, vAt4New] = rv02rvf(rA0New, vA0New, ...
                             X(5), muSunNew);               % RV of Psyche (t=t4)
[rMt5New, vMt5New] = rv02rvf(rM0New, vM0New, ...
                             X(6), muSunNew);               % RV of Mars (t=t5)
[vt4New, vt51New] = LambSol(rAt4New, rMt5New, ...
                           tAENew, muSunNew);               % Lambert problem 4: A->M

% GA-Transfer: M-->E (t5-t6)
tEM2New = X(7) - X(6);
[rEt6New, vEt6New] = rv02rvf(rE0New, vE0New, ...
                             X(7), muSunNew);               % RV of Mars (t=t6)
[vt53New, vt61New] = LambSol(rMt5New, rEt6New, ...
                             tEM2New, muSunNew);            % Lambert problem 5: M->M

% Arrival: E-->E (t6-t7)
tMENew = X(8) - X(7);
[rEt7New, vEt7New] = rv02rvf(rE0New, vE0New, ...
                             X(8), muSunNew);               % RV of Earth (t=t7)
[vt63New, vt7New] = LambSol(rEt6New, rEt7New, ...
                            tMENew, muSunNew);              % Lambert problem 6: M->E

% GA-3: SOI (t5)
[vt521New, vt522New, ~] = SOI_opt(vt11New, vt13New, vMt5New, muMarsNew, X(9), X(13));

% GA-4: SOI (t6)
[vt621New, vt622New, ~] = SOI_opt(vt21New, vt23New, vEt6New, muEarthNew, X(10), X(14));

% Velocity change
% Impulse 1: t0
dvt0New = vt0New - vEt0New;
dvt0NormNew = norm(dvt0New);
dvt0Norm = dvt0NormNew / vUnit;
if dvt0Norm < 4
    dvt0New = zeros(3, 1);
    dvt0Norm = 0;
    dvt0NormNew = 0;
else
    dvt0Vector = dvt0New / dvt0NormNew;
    dvt0New = dvt0New - 4 * dvt0Vector;
    dvt0Norm = dvt0Norm - 4;
    dvt0NormNew = dvt0Norm * vUnit;
end

% Impulse 2-1 2-2: t1
dvt11New = vt121New - vt11New;
dvt11NormNew = norm(dvt11New);
dvt12New = vt13New - vt122New;
dvt12NormNew = norm(dvt12New);

% Impulse 3-1 3-2: t2
dvt21New = vt221New - vt21New;
dvt21NormNew = norm(dvt21New);
dvt22New = vt23New - vt222New;
dvt22NormNew = norm(dvt22New);

% Impulse 4: t3
dvt3New = vAt3New - vt3New;
dvt3NormNew = norm(dvt3New);

% Impulse 5: t4
dvt4New = vt4New - vAt4New;
dvt4NormNew = norm(dvt4New);

% Impulse 6-1 6-2: t5
dvt51New = vt521New - vt51New;
dvt51NormNew = norm(dvt51New);
dvt52New = vt53New - vt522New;
dvt52NormNew = norm(dvt52New);

% Impulse 7-1 7-2: t6
dvt61New = vt621New - vt61New;
dvt61NormNew = norm(dvt61New);
dvt62New = vt63New - vt622New;
dvt62NormNew = norm(dvt62New);

% Impulse 8: t7
dvt7New = vEt7New - vt7New;
dvt7NormNew = norm(dvt7New);
dvt7Norm = dvt7NormNew / vUnit;
if dvt7Norm < 4
    dvt7New = zeros(3, 1);
    dvt7Norm = 0;
    dvt7NormNew = 0;
else
    dvt7Vector = dvt7New / dvt7NormNew;
    dvt7New = dvt7New - 4 * dvt7Vector;
    dvt7Norm = dvt7Norm - 4;
    dvt7NormNew = dvt7Norm * vUnit;
end

% First phase end, total mass loss during 1st phase
dvBeforeSampNormNew = dvt0NormNew + dvt11NormNew + dvt12NormNew + dvt21NormNew + dvt22NormNew + dvt3NormNew;
[mTotalBeforeSamp, dmBeforeSamp] = impulseFuel(mTotal0, dvBeforeSampNormNew, IspNew, g0New);

% Total mass loss after sampling
dvAfterSampNew = dvt4NormNew + dvt51NormNew + dvt52NormNew + dvt61NormNew + dvt62NormNew + dvt7NormNew;
dmAfterSamp = mFuel - dmBeforeSamp;
[~, temp] = impulseFuel(1 / dmAfterSamp, dvAfterSampNew, IspNew, g0New);
mSample = 1 / temp - mTotalBeforeSamp;

J = -mSample;

end