%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOI: Two impulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vMid1, vMid2, dvGA] = SOI_opt(v1, v2, vPlanet, muPlanet, rp, phi)
vMid0 = v1;                                                                 % Use v1 as initial value
obj_func = @(vMid)SOI_impulse(vMid, v1, v2, vPlanet, muPlanet, rp, phi);
options = optimoptions("fminunc", "Display", "off");
[vMid1, ~] = fminunc(obj_func, vMid0, options);
[vMid2, dvGA] = SOI(vMid1, vPlanet, muPlanet, rp, phi);
end

% vMid: velocity right before GA
function dv = SOI_impulse(vMid, v1, v2, vPlanet, muPlanet, rp, phi)
[vMid_after, ~] = SOI(vMid, vPlanet, muPlanet, rp, phi);
dv1 = vMid - v1;
dv2 = v2 - vMid_after;
dv = norm(dv1) + norm(dv2);
end