% Code to generate Figure 4 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagpts_new, mycolours)

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

%% Parameters
nvals = logspace(1,3)'; 
N = numel(nvals);

largestL0  = zeros(N,1);
largestL1  = zeros(N,1);
smallestL0 = zeros(N,1);
smallestL1 = zeros(N,1);

%% Compute roots
for k = 1:N
    n = nvals(k);

    % Golub-Welsch routine to compute roots
    rootsL0 = lagpts_new(n, [0, inf], 'gw', 0);
    rootsL1 = lagpts_new(n, [0, inf], 'gw', 1);

    largestL0(k)  = max(rootsL0);
    largestL1(k)  = max(rootsL1);
    smallestL0(k) = min(rootsL0);
    smallestL1(k) = min(rootsL1);
end


%% ---- Smallest zero ----
figure(1)
loglog(nvals, smallestL0, 's', 'MarkerSize', 8); hold on
loglog(nvals, smallestL1, 'o', 'MarkerSize', 4);
loglog(nvals, 1./nvals, 'k--'); hold off

grid on
ylim([5e-4, 1]);

legend( ...
    '$L_n(x)$', ...
    '$L^{(1)}_n(x)$', ...
    'Location','northeast', 'FontSize', 16);

title('Smallest zero', 'FontSize', 20);
xlabel('$n$', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig4a.pdf', 'ContentType', 'vector');
cd ../figures

%% ---- Largest zero ----
figure(2)
loglog(nvals, largestL0, 's', 'MarkerSize', 8); hold on
loglog(nvals, largestL1, 'o', 'MarkerSize', 4);
loglog(nvals, 4*nvals, 'k--'); hold off

grid on

legend( ...
    '$L_n(x)$', ...
    '$L^{(1)}_n(x)$', ...
    'Location','northwest', 'FontSize', 16);

title('Largest zero', 'FontSize', 20);
xlabel('$n$', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig4b.pdf', 'ContentType', 'vector');
cd ../figures