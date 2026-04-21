function [x, D, w, v] = lagdif_stable(N, M, b)
% LAGDIF_STABLE - Compute Gauss-Laguerre points and differentiation matrix
%   X = LAGDIF_STABLE(N) returns the N-vector [0 ; XL], where XL are the N-1
%   Gauss-Laguerre points.
%
%   [X, D] = LAGDIF_STABLE(N) returns the corresponding first order 
%   differentiation matrix.
%
%   [X, D] = LAGDIF_STABLE(N, 2) returns the corresponding 2nd order 
%   differentiation matrix.
%
%   [X, D, WL, VL] = LAGDIF_STABLE(N) returns also the corresponding barycentric
%   weights, VL, and quadrature weights, W. (Note W are not precisely the
%   optimal quadrature weights for these nodes; they are [0, WL], where
%   WL are the Gauss-Laguerre quadrature weights corresponding to XL.)
%
%   [X, D, W, V] = LAGDIF_STABLE(N, M, B) incorporates the scaling X->X/B in all
%   the outputs.

% Emma Nel, Nick Hale, Stellenbosch University, 2026

if ( nargin == 1 )
    M = 1;
end

if M > 2
    error('Error: Higher order differentiation matrices not implemented.')
end

[x, w, ~, ders] = lagpts_new(N-1, 'glr');                      
v = [1; exp(-x/2)./(x.*ders)];
x = [0; x]; w = [0, w];                                     % Add pt at 0

% First-order matrix
c = [1; ders.*x(2:end)];
C = c./c';                
ii = 1:N+1:N^2;                       
X = 1./(x-x');              
X(ii) = -.5;
D = X.*C;                                                   % Off-diagonals
D(ii) = [0.5-N;0.5./x(2:end)];                              % Diagonals (direct formula)

% Second-order matrix
if M == 2
    D = 2*X.*(C.*repmat(diag(D),1,N) - D);                  % Off-diagonals
    D(ii) = [(2*N^2-2*N+1)/4;(x(2:end).^2 ...               % Diagonals (direct formula)
        + 2*x(2:end) - 4*N*x(2:end) - 4)./(12*x(2:end).^2)];
end

if ( nargin > 2 )                                           % Scaling
    x = x/b; D = (b^M)*D; w = w/b;
end

end