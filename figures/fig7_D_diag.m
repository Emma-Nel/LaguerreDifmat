% Code to generate Figure 7 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagpts_new, mycolours)
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
N = numel(nvals);

% Preallocate
err_DIRECT_glr = zeros(1,N);
err_DMS_dms    = zeros(1,N);
err_DMS_glr    = zeros(1,N);

%% Main loop
for k = 1:N
    n = nvals(k);

    ii = 1:n+2:(n+1)^2; % diagonal indices

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))
    D_diag_hp = D_hp(ii);

    % ----- direct formula, glr roots -----
    x_glr = lagpts_new(n, 'glr');
    D_diag_DIRECT_glr = [0.5-(n+1); 0.5./x_glr];

    % ----- DMSUITE method, dms roots -----
    [~, D] = lagdif(n+1,1,1);
    D_diag_DMS_dms = D(ii)';

    % ----- DMSUITE method, glr roots -----
    x_glr = [0; x_glr]; 
    alpha = exp(-x_glr./2); beta = -0.5*ones(1,n+1);
    D = poldif(x_glr,alpha,beta);
    D_diag_DMS_glr = D(ii)';

    % error
    err_DIRECT_glr(k) = norm(abs(D_diag_DIRECT_glr - D_diag_hp) ./ abs(D_diag_hp), inf);
    err_DMS_dms(k)    = norm(abs(D_diag_DMS_dms - D_diag_hp) ./ abs(D_diag_hp), inf);
    err_DMS_glr(k)    = norm(abs(D_diag_DMS_glr - D_diag_hp) ./ abs(D_diag_hp), inf);
end

% error profile
errN_DIRECT = abs(D_diag_DIRECT_glr - D_diag_hp) ./ abs(D_diag_hp);
errN_DMS    = abs(D_diag_DMS_glr - D_diag_hp) ./ abs(D_diag_hp);

%% -------- Plot: error vs n --------
figure(1)
C = mycolours();
semilogy(nvals, err_DMS_dms, '-o', 'Color', C(1,:)); hold on
semilogy(nvals, err_DMS_glr, '-', 'Color', C(1,:));
semilogy(nvals, err_DIRECT_glr, '-o', 'Color', C(2,:)); hold off
grid on

legend('DMSUITE','DMSUITE (glr roots)', 'direct', ...
    'Location','northwest','FontSize',16)

xlabel('$n$', 'FontSize', 18);
ylabel('inf norm error', 'FontSize', 18);

ylim([1e-16,1e-8])

cd ../output
exportgraphics(gcf, 'fig7a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: error profile n = 500 --------
figure(2)
C = mycolours();
semilogy(1:n+1, errN_DMS, '-','Color', C(1,:)); hold on
semilogy(1:n+1, errN_DIRECT, '-o','Color', C(2,:)); hold off
grid on

legend('DMSUITE (glr roots)','direct', ...
    'Location','northwest','FontSize',16)

xlabel('index $k$ for entry $D_{k,k}$', 'FontSize', 18);

ylim([1e-16,1e-8])
xlim([0,500])

cd ../output
exportgraphics(gcf, 'fig7b.pdf', 'ContentType', 'vector');
cd ../figures