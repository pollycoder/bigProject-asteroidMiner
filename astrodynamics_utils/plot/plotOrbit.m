%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the orbit from orbit elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotOrbit(coe0, mu)
n=100;
tol = 1e-20;
f = linspace(0, 2 * pi, n);
for i = 1:length(f)
    coe = coe0;
    coe(6) = f(i);
    r(:, i) = coe2rv(coe, mu, tol);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', 1.5, 'Color', 'k');hold on
end