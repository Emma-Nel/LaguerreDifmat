% Code to generate Figure 8 of
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
N = numel(nvals);

% Preallocate
err_DMS = zeros(1,N);
err_REC = zeros(1,N);
err_HY  = zeros(1,N);
err_GLR = zeros(1,N);

%% Main loop (first-order)
for k = 1:N
    n = nvals(k);

    ii = 1:n+2:(n+1)^2; % diagonal indices

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))
    D_hp(ii) = [];

    % ----- DMSUITE method -----
    [~, D_dms] = lagdif(n+1, 1, 1);
    D_dms = D_dms(:); D_dms(ii) = [];

    % ----- GLR method -----
    [x_glr, D_glr] = lagdif_stable(n+1, 1, 1);
    D_glr = D_glr(:); D_glr(ii) = [];

    % ----- Recurrence -----
    x_glr(1) = [];
    [~,ders_rec] = laguerre_fun(n,0,x_glr);
    D_rec = lagdif_stable_extract(x_glr, ders_rec);
    D_rec = D_rec(:); D_rec(ii) = [];

    % ----- Huang & Yu -----
    [~,ders_hy] = laguerre_fun_modified_adaptive(n+1,0,x_glr);
    D_hy = lagdif_stable_extract(x_glr,ders_hy);
    D_hy = D_hy(:); D_hy(ii) = [];

    % error
    err_DMS(k) = norm(abs(D_dms - D_hp) ./ abs(D_hp), inf);
    err_REC(k) = norm(abs(D_rec - D_hp) ./ abs(D_hp), inf);
    err_HY(k) = norm(abs(D_hy - D_hp) ./ abs(D_hp), inf);
    err_GLR(k) = norm(abs(D_glr - D_hp) ./ abs(D_hp), inf);

end

%% -------- Plot: Off-diagonals (first-order matrix) --------
err_REC(38) = NaN; % remove blow-up in final finite entry

figure(1)
C = mycolours();
semilogy(nvals, err_DMS, 'o-', 'Color', C(1,:), 'MarkerSize', 6); hold on
semilogy(nvals, err_REC, 'o-', 'color', C(3,:) ,'MarkerSize', 6);
semilogy(nvals, err_HY, 'o-', 'color', C(4,:) ,'MarkerSize', 6);
semilogy(nvals, err_GLR, 'o-', 'color', C(2,:) ,'MarkerSize', 6); hold off
grid on

legend('DMSUITE','Recurrence','Huang and Yu', 'GLR', ...
    'Location','northwest','FontSize',16)

xlabel('$n$', 'FontSize', 18)
ylabel('inf norm error','FontSize', 18)

ylim([1e-16,1e-8])

cd ../output
exportgraphics(gcf, 'fig8a.pdf', 'ContentType', 'vector');
cd ../figures

%% Main loop (second-order)

% Preallocate
err_DMS_dms    = zeros(1,N);
err_DMS_glr    = zeros(1,N);
err_GLR        = zeros(1,N);

for k = 1:N
    n = nvals(k);

    % diagonal indices
    ii = 1:n+2:(n+1)^2;

    % ----- Load high-precision data -----
    load(sprintf('hp_data_n=%d.mat',n))
    D2_hp(ii) = [];

    % ----- GLR method -----
    [x_glr,D2_glr] = lagdif_stable(n+1,2,1);
    D2_glr(ii) = []; D2_glr = D2_glr';

    % ----- DMSUITE (dms roots) -----
    [~,D] = lagdif(n+1,2,1);
    D2_DMS_dms = D(:,:,2);
    D2_DMS_dms(ii) = []; D2_DMS_dms = D2_DMS_dms';

    % ----- DMSUITE (glr roots) -----
    alpha = exp(-x_glr./2);
    beta = [-0.5; (-0.5)^2] .* ones(2,n+1);
    D = poldif(x_glr,alpha,beta);
    D2_DMS_glr = D(:,:,2);
    D2_DMS_glr(ii) = []; D2_DMS_glr = D2_DMS_glr';

    % error
    err_GLR(k) = norm(abs(D2_glr - D2_hp) ./ abs(D2_hp),inf);
    err_DMS_dms(k) = norm(abs(D2_DMS_dms - D2_hp) ./ abs(D2_hp),inf);
    err_DMS_glr(k) = norm(abs(D2_DMS_glr - D2_hp) ./ abs(D2_hp),inf);

end
%% -------- Plot: Off-diagonals (second-order matrix) --------
figure(2)
C = mycolours();
semilogy(nvals, err_DMS_dms, 'o-', 'Color', C(1,:), 'MarkerSize', 6); hold on
semilogy(nvals, err_DMS_glr, '-', 'Color', C(1,:), 'MarkerSize', 6);
semilogy(nvals, err_GLR, 'o-', 'color', C(2,:) ,'MarkerSize', 6); hold off
grid on

legend('DMSUITE','DMSUITE (glr roots)','GLR', 'Location','northwest','FontSize',16)

xlabel('$n$','FontSize', 18)
ylabel('inf norm error','FontSize', 18)

ylim([1e-16,1e-8])

cd ../output
exportgraphics(gcf, 'fig8b.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Functions --------
function D = lagdif_stable_extract(x, ders)
% Extract from LAGDIF_STABLE; uses precomputed nodes and derivatives
x = [0;x];
N = length(x);

c = [1; ders.*x(2:end)];
C = c./c';
ii = 1:N+1:N^2;
X = 1./(x-x');
X(ii) = -.5;
D = X.*C;                             % Off-diagonal entries
D(ii) = [0.5-N;0.5./x(2:end)];        % Diagonal entries

end

function [L_hat,L_hat_deriv] = laguerre_fun(n,alpha,x)
% Standard recurrence relation for generalized Laguerre functions and their
% derivatives evaluated at x
if n == 0
    L_hat = exp(-x/2);
    L_hat_deriv = -0.5*L_hat;
    return;
elseif n == 1
    L_hat = (1 + alpha - x).*exp(-x/2);
    L_hat_deriv = -0.5.*L_hat - exp(-x/2).*ones(size(x));
    return;
end

L_nm2 = exp(-x/2);                                           % L_0^alpha(x)
L_nm1 = (1 + alpha - x).*exp(-x/2);                          % L_1^alpha(x)
L_nm1_deriv = -0.5.*L_nm1 - exp(-x/2);                       % L'_1^alpha(x)

for n = 1:n-1
    L_n = ((2*n + 1 + alpha - x) .* L_nm1 ...               % function
        - (n + alpha) .* L_nm2) ./ (n + 1);
    L_n_deriv = L_nm1_deriv - 0.5*L_nm1 - 0.5*L_n;          % derivative

    % save
    L_nm1_deriv = L_n_deriv;
    L_nm2 = L_nm1;
    L_nm1 = L_n;
end

L_hat = L_nm1;
L_hat_deriv = L_nm1_deriv;

end

function [L_hat,L_hat_deriv] = laguerre_fun_modified_adaptive(n,alpha,x)
% Adaptive procedure introduced by [1] that employs the MODIFIED recurrence 
% relation to evaluate the generalized Laguerre function and derivatives at x
%
% [1] S. Huang and H. Yu, "Improved Laguerre spectral methods with less round-off
% errors and better stability", https://arxiv.org/abs/2212.13255

K1 = 32;                                                 % heuristic values
K2 = 32;

if n == 0
    L_hat = exp(-x/2);
    L_hat_deriv = -0.5*L_hat;
    return;
elseif n == 1
    L_hat = (1 + alpha - x).*exp(-x/2);
    L_hat_deriv = -0.5.*L_hat - exp(-x/2).*ones(size(x));
    return;
end

xb = x/2;
L_n = 1 + alpha - x;                                        % L_1^alpha(x)
delL_n = alpha-x;
L_n_deriv = -ones(size(x));                                 % L'_1^alpha(x)

for n = 1:n-1
    % update
    L_n_deriv = L_n_deriv - L_n;                            % derivative
    delL_n = (n+alpha)./(n+1)*delL_n - x./(n+1).*L_n;
    L_n = L_n+delL_n;                                       % function

    % adaptive weighting
    i = abs(L_n)>exp(K1);
    xc = zeros(size(x));
    if n == 1
        i(:) = true;
    end
    xc(i) = min(max(log(abs(L_n(i)))+K2,0),xb(i));
    delL_n = delL_n.*exp(-xc);
    L_n = L_n.*exp(-xc);
    L_n_deriv = L_n_deriv.*exp(-xc);
    xb = xb-xc;
end

L_hat = L_n.*exp(-xb);
L_hat_deriv = L_n_deriv.*exp(-xb);

end