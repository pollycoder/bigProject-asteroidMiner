%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the trajectory from r0, v0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTrajectory(r0, v0, dt, mu, style)
n1=1000;
t = linspace(0, dt, n1);
for i=1:length(t)
    [r(:, i), ~] = rv02rvf(r0, v0, t(i), mu);
end
plot3(r(1,:), r(2,:), r(3,:), 'LineWidth', style.LineWidth, ...
    'Color', style.LineColor, 'LineStyle',style.LineStyle);hold on
plot3(r(1, end), r(2,end), r(3,end), style.pointStyle, 'LineWidth', 2);hold on
text(r(1, end), r(2,end), r(3,end), style.pointText);
end