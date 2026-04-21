% Code to generate Figure 2 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (mycolours)

clear; clc; close all;

%% Paths
addpath('../src');

%% Defaults
set(groot, 'defaultAxesColorOrder', mycolours());
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%% Grid
x = linspace(0, 10, 500);

%% -------- Plot: various alpha --------
n = 4;

figure(1)
for alpha = 3:-1:0
    L = laguerre_poly(n, alpha, x);
    plot(x, L); hold on
end
hold off
grid on

xlim([0, 10]);
ylim([-10, 10]);

legend( ...
    '$L^{(3)}_{4}$', ...
    '$L^{(2)}_{4}$', ...
    '$L^{(1)}_{4}$', ...
    '$L^{(0)}_{4}$', ...
    'Location','southwest', 'FontSize', 16);

ylabel('$L^{(\alpha)}_n(x)$', 'FontSize', 18);
xlabel('$x$', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig2a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: various n --------
figure(2)
for n = 1:4
    L = laguerre_poly(n, 0, x);
    plot(x, L); hold on
end
hold off
grid on

xlim([0, 10]);
ylim([-10, 10]);

legend( ...
    '$L_{1}$', '$L_{2}$', '$L_{3}$', '$L_{4}$', ...
    'Location','southwest', 'FontSize', 16);

ylabel('$L_n(x)$', 'FontSize', 18);
xlabel('$x$', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig2b.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Functions --------
function [L, L_deriv] = laguerre_poly(N, alpha, x)
% Standard recurrence relation for generalized Laguerre polynomials and
% derivatives evaluated at x
if N == 0
    L = ones(size(x));
    L_deriv = zeros(size(x));
    return;
elseif N == 1
    L = 1 + alpha - x;
    L_deriv = -ones(size(x));
    return;
end

L_nm2 = ones(size(x));              % L_0
L_nm1 = 1 + alpha - x;              % L_1
L_nm1_deriv = -ones(size(x));

for k = 1:N-1
    L_n = ((2*k + 1 + alpha - x).*L_nm1 - (k + alpha).*L_nm2) / (k + 1);
    L_n_deriv = L_nm1_deriv - L_nm1;

    L_nm2 = L_nm1;
    L_nm1 = L_n;
    L_nm1_deriv = L_n_deriv;
end

L = L_nm1;
L_deriv = L_nm1_deriv;

end