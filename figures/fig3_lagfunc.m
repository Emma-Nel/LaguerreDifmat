% Code to generate Figure 3 of
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
x = linspace(0, 100, 1001);

%% -------- Plot: Laguerre polynomials --------
figure(1); hold on
for n = 10:10:30
    L = laguerre_poly(n, 0, x);
    plot(x, L);
end
hold off
grid on

legend({'$L_{10}$', '$L_{20}$', '$L_{30}$'}, ...
    'Location','southwest', 'FontSize', 16);

xlabel('$x$', 'FontSize', 18);
ylabel('$L_n(x)$', 'FontSize', 18);
ylim([-1e10, 1e10]);

cd ../output
exportgraphics(gcf, 'fig3a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: Laguerre functions --------
figure(2); hold on
for n = 10:10:30
    L_hat = exp(-x/2) .* laguerre_poly(n, 0, x);
    plot(x, L_hat);
end
hold off
grid on

legend({'$\widehat{L}_{10}$', '$\widehat{L}_{20}$', '$\widehat{L}_{30}$'}, ...
    'Location','southwest', 'FontSize', 16);

xlabel('$x$', 'FontSize', 18);
ylabel('$\widehat{L}_n(x)$', 'FontSize', 18);
ylim([-1, 1]);

cd ../output
exportgraphics(gcf, 'fig3b.pdf', 'ContentType', 'vector');
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