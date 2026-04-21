% Code to generate Figure 1 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagdif_stable, lagpts_new, mycolours)
%   - data (high-precision data files)
%   - dmsuite (lagdif, poldif, lagroots)

clear; clc; close all;

%% Paths
addpath('../src');
addpath('../data');

%% Check dependencies

files = {'lagdif.m','poldif.m','lagroots.m'};

if any(~cellfun(@(f) exist(f,'file'), files))
    error('Error: This code requires <a href="https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite">DMSUITE</a> to run.')
end

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
errinf_GLR_offdiag = zeros(1, N);
errinf_GLR_diag    = zeros(1, N);
errinf_DMS_offdiag = zeros(1, N);
errinf_DMS_diag    = zeros(1, N);

%% Main loop
for k = 1:N
    n = nvals(k);

    ii = 1:n+2:(n+1)^2; % diagonal indices

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))

    % Split diagonal/off-diagonal
    D_hp_diag = D_hp(ii); D_hp(ii)  = [];

    % ----- GLR method -----
    [~, D_glr] = lagdif_stable(n+1, 1, 1);
    D_glr = D_glr(:);
    D_glr_diag = D_glr(ii); D_glr(ii)  = [];

    errinf_GLR_offdiag(k) = norm(abs(D_glr - D_hp) ./ abs(D_hp), inf);
    errinf_GLR_diag(k)    = norm(abs(D_glr_diag - D_hp_diag) ./ abs(D_hp_diag), inf);

    % ----- DMSUITE method -----
    [~, D_dms] = lagdif(n+1, 1, 1);
    D_dms = D_dms(:);
    D_dms_diag = D_dms(ii); D_dms(ii)  = [];

    errinf_DMS_offdiag(k) = norm(abs(D_dms - D_hp) ./ abs(D_hp), inf);
    errinf_DMS_diag(k)    = norm(abs(D_dms_diag - D_hp_diag) ./ abs(D_hp_diag), inf);
end

%% -------- Plot: Diagonals --------
figure(1)
semilogy(nvals, errinf_DMS_diag, 'o-', 'MarkerSize', 6); hold on
semilogy(nvals, errinf_GLR_diag, 'o-', 'MarkerSize', 6); hold off
grid on

ylim([1e-16,1e-8])

legend('DMSUITE', 'Proposed method', ...
    'FontSize', 16, 'Location', 'northwest');

xlabel('$n$','FontSize', 18);
ylabel('inf norm error', 'FontSize', 18);
title('Diagonal entries', 'FontSize', 20);

cd ../output
exportgraphics(gcf, 'fig1a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: Off-diagonals --------
figure(2)
semilogy(nvals, errinf_DMS_offdiag, 'o-', 'MarkerSize', 6); hold on
semilogy(nvals, errinf_GLR_offdiag, 'o-', 'MarkerSize', 6); hold off
grid on

ylim([1e-16,1e-8])

legend('DMSUITE', 'Proposed method', ...
    'FontSize', 16, 'Location', 'northwest');

xlabel('$n$','FontSize', 18);
ylabel('inf norm error', 'FontSize', 18);
title('Off-diagonal entries', 'FontSize', 20);

cd ../output
exportgraphics(gcf, 'fig1b.pdf', 'ContentType', 'vector');
cd ../figures