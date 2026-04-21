% Code to generate Figure 9 of
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
N    = numel(nvals);

% Preallocate
err_diag_DIRECT_glr = zeros(1,N);
err_diag_DIRECT_hp  = zeros(1,N);
err_diag_DMS_dms    = zeros(1,N);
err_diag_DMS_glr    = zeros(1,N);

%% Main loop
for k = 1:N
    n = nvals(k);

    % diagonal indices
    ii = 1:n+2:(n+1)^2;

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))

    % split diagonal/off-diagonal
    D2_diag = D2_hp(ii);

    % ----- direct (glr roots) -----
    x_glr = lagpts_new(n, 'glr');
    d_diag_DIRECT_glr = [(2*(n+1)^2-2*(n+1)+1)/4; ...
        (2*x_glr - 4*(n+1)*x_glr + x_glr.^2 - 4)./(12*x_glr.^2)];

    % ----- direct (hp roots) -----
    d_diag_DIRECT_hp = [(2*(n+1)^2-2*(n+1)+1)/4; ...
        (2*x_hp - 4*(n+1)*x_hp + x_hp.^2 - 4)./(12*x_hp.^2)];

    % ----- DMSUITE (dms roots) -----
    [~,D] = lagdif(n+1,2,1);
    D2 = D(:,:,2);
    d_diag_DMS_dms = D2(ii)';

    % ----- DMSUITE (glr roots) -----
    x_glr = [0; x_glr];
    alpha = exp(-x_glr./2);
    beta = [-0.5; (-0.5)^2] .* ones(2,n+1);
    D = poldif(x_glr,alpha,beta);
    d_diag_DMS_glr = D(:,:,2);
    d_diag_DMS_glr = d_diag_DMS_glr(ii)';

    % error
    err_diag_DIRECT_glr(k) = norm(abs(d_diag_DIRECT_glr - D2_diag)./abs(D2_diag), inf);
    err_diag_DIRECT_hp(k)  = norm(abs(d_diag_DIRECT_hp - D2_diag)./abs(D2_diag), inf);
    err_diag_DMS_dms(k)    = norm(abs(d_diag_DMS_dms - D2_diag)./abs(D2_diag), inf);
    err_diag_DMS_glr(k)    = norm(abs(d_diag_DMS_glr - D2_diag)./abs(D2_diag), inf);
end

% error profile
errN_dms_rel    = (abs(d_diag_DMS_glr - D2_diag) ./ abs(D2_diag));
errN_direct_rel = (abs(d_diag_DIRECT_glr - D2_diag) ./ abs(D2_diag));
errN_hp_rel     = (abs(d_diag_DIRECT_hp - D2_diag) ./ abs(D2_diag));

%% -------- Plot: error vs n --------
figure(1)
C = mycolours();
semilogy(nvals, err_diag_DMS_dms, 'o-', 'color', C(1,:), 'MarkerSize', 6); hold on
semilogy(nvals, err_diag_DMS_glr, '-', 'color', C(1,:),'MarkerSize', 6);
semilogy(nvals, err_diag_DIRECT_glr, 'o-', 'color', C(2,:),'MarkerSize', 6); hold off
grid on

ylim([1e-16,1e-8])

legend('DMSUITE', 'DMSUITE (glr roots)', ...
    'Direct','FontSize',16,'Location','northwest');

xlabel('$n$', 'FontSize', 18);
ylabel('inf norm error', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig9a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Plot: error profile n = 500 --------
figure(2)
semilogy(1:n+1, errN_dms_rel, '-', 'MarkerSize', 6); hold on
semilogy(1:n+1, errN_direct_rel, '-o', 'MarkerSize', 5);
semilogy(1:n+1, errN_hp_rel, '-o', 'MarkerSize', 5); hold off
grid on

xlim([0,500])
ylim([1e-16,1e-8])

legend('DMSUITE (glr roots)', 'direct', 'direct (hp roots)','FontSize',16,'Location','northwest');

xlabel('index $k$ for entry $D_{k,k}$', 'FontSize', 18);

cd ../output
exportgraphics(gcf, 'fig9b.pdf', 'ContentType', 'vector');
cd ../figures