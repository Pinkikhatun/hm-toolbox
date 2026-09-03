%TEST_HSS_EUCLIDEAN_GRADIENT
clear; 
% rng(0);

r = 2;
hssoption('block-size', 4);

n = 16;                            % gives a 3-level tree with block-size 4
H = hssgallery('rand', n, r);      % current iterate X
A = randn(n);                      % target matrix

f = @(K) 0.5 * norm(full(K) - A, 'fro')^2;

%% ---- 1. gradient -------------------------------------------------------
G = hss_euclidean_gradient_mf(H, A);

%% 
TH = hsstangent('randn', H);       % arbitrary direction in generator space

 ana = hsstangent_inner(G, TH);     % <grad f, TH>
% ana = dot(G, TH);

fprintf('\n   t          FD quotient        ana         rel.err\n');
fprintf('  ------------------------------------------------------------\n');
for t = 10.^-(3:8)
    num = ( f(hss_perturb(H, TH,  t)) - f(hss_perturb(H, TH, -t)) ) / (2*t);
    fprintf('  %8.1e   %+14.8e   %+14.8e   %8.2e\n', ...
        t, num, ana, abs(num - ana) / max(1, abs(ana)));
end

%% 
%  <grad_C , C> summed over all blocks + <grad_D, D> = <E, X>
E = full(H) - A;
lhs = euler_check(H, G);
rhs = E(:).' * reshape(full(H), [], 1);
fprintf('\n  Euler identity  |<grad,theta>_mult - <E,X>| = %8.2e\n', ...
    abs(lhs - rhs));

% (b) gradient of an exact representation is zero
He = hss(full(H));                 % same matrix, re-compressed
Ge = hss_euclidean_gradient_mf(He, full(He));
fprintf('  ||grad f(H,H)||                                    = %8.2e\n', ...
    sqrt(hsstangent_inner(Ge, Ge)));

%% ---- 4. Riemannian gradient ---------------------------------
 rgrad = projection(G);          
% -------------------------------------------------------------------------
function s = euler_check(H, G)
%EULER_CHECK  
if H.leafnode
    s = G.TD(:).' * H.D(:);
else
    dU = H.U;  %
    s = G.TB12(:).' * H.B12(:) + G.TB21(:).' * H.B21(:) ...
        + euler_check(H.A11, G.TA11) + euler_check(H.A22, G.TA22);
end
end
