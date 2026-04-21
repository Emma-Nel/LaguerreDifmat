% Code to generate Figure 6 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagpts_new, mycolours)
%   - data (high-precision data files)
%   - dmsuite (lagroots)

clear; clc; close all;

%% Paths
addpath('../src');
addpath('../data');

%% Check dependencies

if ( ~exist('lagroots.m','file') ), error(['Error: This code requires ' ...
        '<a href="https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite">DMSUITE</a> to run.']), end

%% Defaults
set(groot, 'defaultAxesColorOrder', mycolours());
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

%% Parameters
nvals = 10:10:500;
N     = numel(nvals);

% Preallocate
err_roots_GLR_rel = zeros(1,N);
err_roots_DMS_rel = zeros(1,N);

%% Main loop
for k = 1:N
    n = nvals(k);

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))

    % ----- GLR method -----
    x_glr = lagpts_new(n, 'glr');  
    err_roots_GLR_rel(k) = norm(abs(x_glr - x_hp) ./ abs(x_hp), inf);

    % ----- DMSUITE method -----
    x_dms = lagroots(n);
    err_roots_DMS_rel(k) = norm(abs(x_dms - x_hp) ./ abs(x_hp), inf);
end

%% -------- Plot: error vs n --------
figure(1)
semilogy(nvals, err_roots_DMS_rel, '-'); hold on
semilogy(nvals, err_roots_GLR_rel, '-'); hold off
grid on

ylim([1e-16,1e-8])

xlabel('$n$', 'FontSize', 18);
ylabel('inf norm error', 'FontSize', 18);

legend('DMSUITE $\,$','GLR $\,$','Location','northwest','FontSize',16)

cd ../output
exportgraphics(gcf, 'fig6a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: error profile n = 500 --------
err_glr     = abs(x_glr - x_hp) ./ abs(x_hp);
err_dmsuite = abs(x_dms - x_hp) ./ abs(x_hp);

figure(2)
semilogy(1:n, err_dmsuite, '.-','MarkerSize',15); hold on
semilogy(1:n, err_glr, '.-','MarkerSize',15); hold off
grid on

ylim([1e-16,1e-8])

xlabel('index $k$ for root $x_k$', 'FontSize', 18);

legend('DMSUITE $\,$','GLR $\,$','Location','northwest','FontSize',16)

cd ../output
exportgraphics(gcf, 'fig6b.pdf', 'ContentType', 'vector');
cd ../figures