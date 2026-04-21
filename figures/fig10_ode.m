% Code to generate Figure 10 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagdif_stable, lagpts_new, mycolours)
%   - dmsuite (lagdif, poldif, lagroots)

clear; clc; close all;

%% Paths
addpath('../src');

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
nvals = 10:2:250;
N = numel(nvals);
gamma = 2;
beta = 4.03;

%% Exact solution 
u_exact = @(x) exp(-x/4) .* sin(2*x);
f = @(x) (exp(-x/4) .* (16*cos(2*x) + 95*sin(2*x))) ./ 16;

lbc = u_exact(0);

% Preallocate 
err_glr = zeros(1,N);
err_dms  = zeros(1,N);

%% Main loop 
for k = 1:N
    n = nvals(k);
    fprintf('N = %d\n', n);

    I = eye(n);

    % ----- Proposed method -----
    [x_glr,D2_glr] = lagdif_stable(n, 2, beta);
    A = -D2_glr + gamma*I;

    rhs = f(x_glr);
    A(1,:) = 0; A(1,1) = 1;                 % row replacement for bc
    rhs(1) = lbc;

    u = A \ rhs;
    err_glr(k) = norm(abs(u - u_exact(x_glr)), inf);

    % ----- DMSUITE -----
    [x_dms,D2_dms] = lagdif(n, 2, beta);
    D2_dms = D2_dms(:,:,2);

    A = -D2_dms + gamma*I;
    rhs = f(x_dms);

    A(1,:) = 0; A(1,1) = 1;                 % row replacement for bc
    rhs(1) = lbc;

    u = A \ rhs;
    err_dms(k) = norm(abs(u - u_exact(x_dms)), inf);
end

%% ----- Plot: solution -----
x_plot = linspace(0,20,1000);
x_nodes = lagpts_new(230, 'glr'); x_nodes = x_nodes./beta;

figure(1)
plot(x_plot, u_exact(x_plot), 'k'); hold on
plot(x_nodes, u_exact(x_nodes), 'o', 'MarkerSize', 4); hold off
grid on

xlim([0,20])
ylim([-1,1])

xlabel('$x$','FontSize', 18)
ylabel('$u(x)$','FontSize', 18)

cd ../output
exportgraphics(gcf, 'fig10a.pdf', 'ContentType', 'vector');
cd ../figures

%% ----- Plot: convergence -----
figure(2)
semilogy(nvals, err_dms, '-'); hold on
semilogy(nvals, err_glr, '--'); hold off
grid on

legend('DMSUITE', 'Proposed method', 'Location', 'southwest','FontSize', 16)

xlabel('$n$','FontSize', 18)
ylabel('inf norm error','FontSize', 18)

cd ../output
exportgraphics(gcf, 'fig10b.pdf', 'ContentType', 'vector');
cd ../figures
