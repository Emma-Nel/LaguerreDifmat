% Code to generate Figure 5 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (mycolours)
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

%% -------- Figure 5 (left): min/max entries of D --------
nvals = 10:150;
N = numel(nvals);

minD = NaN(1,N);
maxD = NaN(1,N);

for k = 1:N
    n = nvals(k);

    [~, D] = lagdif(n,1,1);

    if any(isnan(D(:))) || any(isinf(D(:)))
        break
    end

    minD(k) = min(abs(D(:)));
    maxD(k) = max(abs(D(:)));
end

figure(1)
C = mycolours();
semilogy(nvals, maxD, '-','Color',C(1,:)); hold on
semilogy(nvals, minD, '-','Color',C(1,:)); hold off

xlim([50,130])
ylim([1e-4,1e3])
grid on

xlabel('$n$', 'FontSize', 18)
ylabel('magnitude', 'FontSize', 18)

% Annotations
text(85, 1.2e-3, 'min', 'FontSize', 16)
text(85, 4e1, 'max', 'FontSize', 16)

cd ../output
exportgraphics(gcf, 'fig5a.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Figure 5 (right): asymptotic scaling --------
nvals = 10:400;
N = numel(nvals);

vn_inv          = zeros(1,N);
vn_scaled_inv   = zeros(1,N);
cn              = zeros(1,N);
cn_scaled       = zeros(1,N);
E               = zeros(1,N);

for k = 1:N
    n = nvals(k);

    x = lagroots(n);

    % Pairwise differences
    X = x' - x;
    X(1:n+1:end) = 1;

    % Barycentric weights (product form)
    vn_inv(k) = prod(X(:,end));

    % Scaled barycentric weights (derivative form)
    [~, dL] = laguerre_poly(n,0,x);
    vn_scaled_inv(k) = abs(dL(end));

    % Scaled quantities
    cn(k)        = abs(vn_inv(k) * exp(-x(end)/2));
    cn_scaled(k) = abs(vn_scaled_inv(k) * exp(-x(end)/2));

    % Exponential term
    E(k) = abs(exp(-x(end)/2));
end

figure(2)
C = mycolours();
semilogy(nvals, vn_inv, '-','Color',C(2,:)); hold on
semilogy(nvals, vn_scaled_inv, '--','Color',C(2,:));
semilogy(nvals, cn, '-','Color',C(1,:));
semilogy(nvals, cn_scaled, '--','Color',C(1,:));
semilogy(nvals, E, '-','Color',C(3,:));

% Machine limits
semilogy(nvals, realmax*ones(size(nvals)), 'k--', 'LineWidth', 2);
semilogy(nvals, realmin*ones(size(nvals)), 'k--', 'LineWidth', 2); hold off

xlim([10,400])
ylim([1e-310,1e310])
grid on

xlabel('$n$', 'FontSize', 18)
ylabel('magnitude', 'FontSize', 18)

yticks([1e-300,1e-200,1e-100,1e0,1e100,1e200,1e300])

% Annotations
text(130, 1e280, '$1/v_n$', 'FontSize', 18)
text(310, 1e230, '$1/\tilde{v}_n$', 'FontSize', 18)
text(130, 1e200, '$c_n$', 'FontSize', 18)
text(320, 1e30, '$\tilde{c}_n$', 'FontSize', 18)
text(130, 1e-70, '$e^{-x_n/2}$', 'FontSize', 18)
text(15, 1e280, '\texttt{realmax}', 'FontSize', 14)
text(15, 1e-280, '\texttt{realmin}', 'FontSize', 14)

cd ../output
exportgraphics(gcf, 'fig5b.pdf', 'ContentType', 'vector');
cd ../figures

%% -------- Functions --------
function [L, L_deriv] = laguerre_poly(n, alpha, x)
% Standard recurrence relation for generalized Laguerre polynomials and
% derivatives evaluated at x
if n == 0
    L = ones(size(x));
    L_deriv = zeros(size(x));
    return;
elseif n == 1
    L = 1 + alpha - x;
    L_deriv = -ones(size(x));
    return;
end

L_nm2 = ones(size(x));              % L_0
L_nm1 = 1 + alpha - x;              % L_1
L_nm1_deriv = -ones(size(x));

for k = 1:n-1
    L_n = ((2*k + 1 + alpha - x).*L_nm1 - (k + alpha).*L_nm2) / (k + 1);
    L_n_deriv = L_nm1_deriv - L_nm1;

    L_nm2 = L_nm1;
    L_nm1 = L_n;
    L_nm1_deriv = L_n_deriv;
end

L = L_nm1;
L_deriv = L_nm1_deriv;

end