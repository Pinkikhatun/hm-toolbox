%% RUN_HSS_RGD_MANOPT
% Riemannian gradient descent for the matrix-nearness problem
%
%       min_{X in M_HSS}  f(X) = 0.5*||X - A||_F^2
%
% using the HSSFIXEDRANKFACTORY manifold and Manopt's STEEPESTDESCENT.

clear;
clc;

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(this_file));
addpath(repo_root);
addpath(genpath(fullfile(repo_root, 'manopt', 'manopt')));

r = 3;
block_size = 32;
n = 512;
hssoption('block-size', block_size);

H = hssgallery('rand', n, r);      % starting point on M_HSS

% Build the same kernel target as in run_hss_rgd.m.
x = sort(rand(n, 1));
A = exp(-abs(x - x.')/0.1);
A = A + 1e-3*rand(n, 1);
normA = norm(A, 'fro');

fprintf('Ordinary numerical rank of A = %d\n', rank(A));
fprintf(' n = %d, HSS rank r = %d, block-size = %d\n', ...
        n, r, block_size);
fprintf(' f(X_0) = %.8e\n', objective(H, A));

%% Describe the problem with the Manopt interface.
problem.M = hssfixedrankfactory(H);
problem.cost = @(X) objective(X, A);

% This is the Euclidean gradient in the ambient space of HSS generators.
% Manopt calls problem.M.egrad2rgrad, which performs the horizontal
% projection implemented by @hsstangent/projection.m.
problem.egrad = @(X) hss_euclidean_gradient_mf(X, A);

metrics.relres = @(X) norm(full(X) - A, 'fro')/normA;

options.maxiter = 1000;
options.tolgradnorm = 1e-4;
options.minstepsize = 1e-14;
options.verbosity = 2;
options.linesearch = @linesearch;
options.ls_contraction_factor = 0.5;
options.ls_suff_decr = 1e-4;
options.ls_optimism = 1.5;
options.ls_max_steps = 50;
options.statsfun = statsfunhelper(metrics);

% Manopt parameterizes its first trial by the norm of the step.  This
% choice makes the first multiplier of -grad equal to alpha0 = 1, as in
% run_hss_rgd.m.
alpha0 = 1;
G0 = problem.M.egrad2rgrad(H, problem.egrad(H));
options.ls_initial_stepsize = alpha0*problem.M.norm(H, G0);

%% Riemannian gradient descent through Manopt.
[X, final_cost, info] = steepestdescent(problem, H, options);

%% Report the same sanity checks as run_hss_rgd.m.
fprintf('\n  rank of the returned HSS matrix : %d (requested %d)\n', ...
        hssrank(X), r);
fprintf('  ||X - A||_F / ||A||_F           : %.6e\n', ...
        norm(full(X) - A, 'fro')/normA);
fprintf('\n  f(X_final)                      : %.6e\n', final_cost);

Acomp = hss(A);
fprintf('\n  for reference, hss(A) has rank  : %d, rel. error %.3e\n', ...
        hssrank(Acomp), norm(full(Acomp) - A, 'fro')/normA);

% First-order consistency of the final Riemannian gradient.
g = problem.M.egrad2rgrad(X, problem.egrad(X));
g2 = problem.M.inner(X, g, g);
fX = problem.cost(X);
fprintf('\n  gradient consistency at the final point:\n');
fprintf('      t        (f(X)-f(R(-t g)))/(t||g||^2)\n');
for t = 10.^-(2:6)
    ratio = (fX - problem.cost(problem.M.retr(X, g, -t)))/(t*g2);
    fprintf('   %8.1e   %12.6f\n', t, ratio);
end
fprintf('   (should approach 1 as t -> 0)\n');

% The standard Manopt history is available as, for example:
%   [info.cost], [info.gradnorm], [info.stepsize], [info.relres]

function value = objective(X, A)
value = 0.5*norm(full(X) - A, 'fro')^2;
end
