% Code to generate Figure 11 of
% "Construction of Laguerre pseudospectral differentiation matrices"
%
% Emma Nel, Stellenbosch University, Apr 2026
%
% Dependencies:
%   - src (lagdif_stable, lagpts_new, mycolours)

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
beta = 10;
nvals = 200:-2:30;
N = numel(nvals);

% indices of eigenvalues/eigenvectors to display
eig_index = [25:-5:5, 1];
nk = numel(eig_index);

% constants
R = 5.08685476;
a = 0.929852862;

% Preallocate
err = zeros(nk, N);
eigenvals = zeros(1,nk);
eigenvectors = zeros(nvals(1),nk);

%% Main loop

% loop over eigenvalues
for k = 1:nk
    eig_ref_index = eig_index(k);

    e_store = zeros(1,N);

    % loop of N values
    for i = 1:N
        n = nvals(i);

        [x,D] = lagdif_stable(n+1,2,beta);
        D2 = D(2:n+1,2:n+1);
        x = x(2:n+1);

        Q = diag(1./(1+exp((x-R)/a)));              % Woods-Saxon potential

        [V,E] = eig(-D2 + eye(size(D2)), Q);
        E = diag(E);

        [E,idx] = sort(abs(E));                     % sort eigenvalues in order of magnitude
        V = V(:,idx);

        if i == 1                                   % take largest N to be the reference solution
            e_ref = E(eig_ref_index);
            e_store(i) = e_ref;
            eigenvals(k) = e_ref;
            eigenvectors(:,k) = V(:,eig_ref_index);
            x_ref = x;
        else
            [~,ind] = min(abs(E - e_ref));
            e_store(i) = E(ind);
        end
    end

    err(k,:) = abs(e_store - e_ref) ./ abs(e_ref);
end

%% ----- Plot: eigenvectors -----
V = eigenvectors;
xx = linspace(0,20,1000)';

figure(1); hold on
for k = 1:size(V,2)
    V(:,k) = V(:,k) / norm(V(:,k));                 % normalize
    V(:,k) = V(:,k) * sign(V(1,k));                 % fix sign

    s = spline(x_ref, V(:,k), xx);                  % interpolate for smoother eigenmodes
    plot(xx, s*20 + eigenvals(k)); hold on          % rescale and shift
end

f = (1 + exp((xx-R)/a));                            % Woods-Saxo potential
plot(xx, f, 'k'); hold off
grid on

xlim([0,15])
ylim([0, eigenvals(1) + 10])

legend_labels = arrayfun(@(k) sprintf('$v_{%d}$',k), eig_index, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northeast','FontSize', 16)

xlabel('$x$','FontSize', 18)
ylabel('eigenmodes $v_k$ shifted by $\lambda_k$','FontSize', 18)

cd ../output
exportgraphics(gcf, 'fig11a.pdf', 'ContentType', 'vector');
cd ../figures

%% ----- Plot: eigenvalue convergence -----
figure(2)
for k = 1:nk
    semilogy(nvals, err(k,:), '-'); hold on
end
hold off

xlim([30,200])
grid on

legend_labels = arrayfun(@(k) sprintf('$\\lambda_{%d}$',k), eig_index, 'UniformOutput', false);
legend(legend_labels, 'Location', 'northeast','FontSize', 16)

xlabel('$n$','FontSize', 18)
ylabel('relative error','FontSize', 18)

cd ../output
exportgraphics(gcf, 'fig11b.pdf', 'ContentType', 'vector');
cd ../figures