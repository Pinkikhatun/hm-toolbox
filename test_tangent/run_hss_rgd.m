%%*****************************************************************
%% RUN_HSS_RGD
%% Riemannian gradient descent for the matrix nearness problem
%%
%%      min_{X in M_HSS}  f(X) = 0.5*||X - A||_F^2
%%


clear; clear functions; clc;
% addpath('/Users/pinki/Downloads/HSS_MATLAB/hm-toolbox')
% rng(0);

r = 2;
hssoption('block-size', 4);
n = 16;

H = hssgallery('rand', n, r);      % starting point on M_HSS
% A = randn(n);                      % target
% A1 = randn(n, 2);
% A2= randn(n,2);
% A=A1*A2';% target
A = zeros(n,n);

% Random points
x = sort(rand(n,1));

% Kernel matrix
for i = 1:n
    for j = 1:n
        A(i,j) = exp(-abs(x(i)-x(j))/0.1);
    end
end

% Add a small random diagonal perturbation
A = A + 1e-3*diag(rand(n,1));

fprintf('Ordinary numerical rank of A = %d\n', rank(A));

f = @(K) 0.5 * norm(full(K) - A, 'fro')^2;

fprintf(' n = %d, HSS rank r = %d, block-size = %d\n', n, r, 4);
fprintf(' f(X_0) = %.8e\n', f(H));

%% ---------------------------------------------------------------
% check_riemannian_gradient(H, A);

%% ---------------------------------------------------------------
%  2. Riemannian gradient descent
%% ---------------------------------------------------------------
opts            = struct();
opts.maxiter    = 1000;
opts.tol        = 1e-4;
opts.linesearch = 'armijo';    % 'armijo' | 'fixed'
opts.alpha0     = 1;
opts.c1         = 1e-4;
opts.beta       = 0.5;
opts.display    = true;
opts.plots      = true;

[X, info] = hss_rgd(H, A, opts);

%% ---------------------------------------------------------------
%  3. Sanity checks and a HONEST reference value
%% ---------------------------------------------------------------
fprintf('\n  rank of the returned HSS matrix : %d (requested %d)\n', ...
        hssrank(X), r);
fprintf('  ||X - A||_F / ||A||_F           : %.6e\n', ...
        norm(full(X) - A,'fro')/norm(A,'fro'));
% -----------------------------------------------------------------
% [lb, blkranks] = hss_bestblock_bound(H, A);

fprintf('\n  f(X_final)                      : %.6e\n', info.cost(end));
% fprintf('  lower bound on min f            : %.6e\n', lb);
% fprintf('  ratio f_final / lower bound     : %.3f\n', info.cost(end)/lb);
% fprintf('  numerical ranks of the off-diagonal blocks of A: %s\n', ...
        % mat2str(blkranks));

Acomp = hss(A);
fprintf('\n  for reference, hss(A) has rank  : %d, rel. error %.3e\n', ...
        hssrank(Acomp), norm(full(Acomp)-A,'fro')/norm(A,'fro'));
% fprintf('  (that is why the old ''hss(A,''''rank'''',r)'' comparison was\n');
% fprintf('   misleading: it was a rank-%d, not a rank-%d, approximation)\n', ...
        % hssrank(Acomp), r);
