function [X, info] = hss_rgd(H, A, opts)
%HSS_RGD  Riemannian gradient descent for  min_{X in M_r} 0.5*||X - A||_F^2
%   on the fixed-rank HSS manifold.
%
%   [X, info] = HSS_RGD(H, A, opts)
%
%   ITERATION
%       E_k    = X_k - A
%       G_k    = grad^E f(X_k)            (hss_euclidean_gradient_mf)
%       g_k    = P_{T_{X_k} M_r}(G_k)     (projection)
%       X_{k+1}= R_{X_k}(-alpha_k g_k)    (retraction)
%
%    <T1,T2> = sum_nodes <dU1,dU2> + <dV1,dV2> + <dR1,dR2>
%                          + <dW1,dW2> + <dB1,dB2> + <dD1,dD2>,
%       dU = U*TU + PU,  dR = R*TR + PR,  ...
%
%   ARMIJO.  With eta = -g the Riemannian sufficient-decrease condition
%       f(R_X(alpha*eta)) <= f(X) + c1*alpha*<grad f(X), eta>_X
%   becomes
%       f(R_X(-alpha*g)) <= f(X) - c1*alpha*||g||_X^2,
%   which is what the backtracking loop below tests.
%
%   OPTS (all optional):
%     .maxiter   maximum iterations              (default 200)
%     .tol       stop when ||grad||_X < tol      (default 1e-8)
%     .linesearch 'armijo' (default) | 'fixed'
%     .alpha0    initial / fixed step            (default 1)
%     .c1        Armijo constant                 (default 1e-4)
%     .beta      backtracking factor             (default 0.5)
%     .alphamin  give up below this step         (default 1e-14)
%     .verbose   print the iteration table       (default true)
%     .plots     draw the three figures          (default true)
%
%   REQUIRES: hss_euclidean_gradient_mf, projection, retraction,
%             hsstangent_inner, hsstangent, base_point.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'maxiter'),    opts.maxiter    = 200;      end
if ~isfield(opts,'tol'),        opts.tol        = 1e-8;     end
if ~isfield(opts,'linesearch'), opts.linesearch = 'armijo'; end
if ~isfield(opts,'alpha0'),     opts.alpha0     = 1;        end
if ~isfield(opts,'c1'),         opts.c1         = 1e-4;     end
if ~isfield(opts,'beta'),       opts.beta       = 0.5;      end
if ~isfield(opts,'alphamin'),   opts.alphamin   = 1e-14;    end
if ~isfield(opts,'display'),    opts.display    = true;     end
if ~isfield(opts,'plots'),      opts.plots      = true;     end
if ~isfield(opts,'ftol'),       opts.ftol       = 1e-12;    end
if ~isfield(opts,'stagmax'),    opts.stagmax    = 10;       end
if ~isfield(opts,'diagnose'),   opts.diagnose   = true;     end

armijo = strcmpi(opts.linesearch,'armijo');

f     = @(K) 0.5 * norm(full(K) - A, 'fro')^2;
normA = norm(A,'fro');

X       = H;
alpha   = opts.alpha0;
bt_prev = 0;
stag    = 0;

cost   = nan(opts.maxiter+1,1);
gnorm  = nan(opts.maxiter+1,1);
steps  = nan(opts.maxiter+1,1);
relres = nan(opts.maxiter+1,1);

if opts.display
    fprintf(['\n  iter        cost           ||grad||_X        alpha    ' ...
             '  rel.res   bt   model\n']);
    fprintf(['  ------------------------------------------------------' ...
             '-------------------------\n']);
end

exitmsg = 'maxiter reached';
k = 0;

for k = 0:opts.maxiter

    fk = f(X);

    % ---- Euclidean gradient, then Riemannian gradient -----------------
    G     = hss_euclidean_gradient_mf(X, A);
    rgrad = projection(G);              % projection(TH): base point is inside

    gn2 = hsstangent_inner(rgrad, rgrad);
    gn  = sqrt(max(gn2,0));

    cost(k+1)   = fk;
    gnorm(k+1)  = gn;
    relres(k+1) = sqrt(2*fk)/normA;     % ||X-A||_F / ||A||_F

    if gn < opts.tol
        steps(k+1) = 0;
        if opts.display
            fprintf('  %4d   %14.8e   %14.8e   %9s  %9.3e  %2s  %6s\n', ...
                    k, fk, gn, '-', relres(k+1), '-', '-');
        end
        exitmsg = sprintf('converged: ||grad||_X = %.3e < tol = %.3e', ...
                          gn, opts.tol);
        break
    end

   %%
    bt    = 0;
    model = NaN;
    if armijo
        if k > 0 && bt_prev == 0
            alpha = min(1, 1.5*alpha);      % gentle growth 
        end
        while true
            Xtry  = retraction(X, rgrad, -alpha);
            ftry  = f(Xtry);
            model = (fk - ftry) / (alpha*gn2);   % ~1 for small alpha if the
                                                 % gradient is correct
            if ftry <= fk - opts.c1*alpha*gn2
                break
            end
            %
            denom = 2*(ftry - fk + gn2*alpha);
            if denom > 0
                anew = gn2*alpha^2/denom;
                anew = min(max(anew, 0.1*alpha), 0.5*alpha);   % safeguard
            else
                anew = opts.beta*alpha;
            end
            alpha = anew;
            bt    = bt + 1;
            if alpha < opts.alphamin
                exitmsg = 'line search failed (step below alphamin)';
                Xtry    = X;
                break
            end
        end
    else
        alpha = opts.alpha0;
        Xtry  = retraction(X, rgrad, -alpha);
        ftry  = f(Xtry);
        model = (fk - ftry) / (alpha*gn2);
    end
    bt_prev = bt;

    steps(k+1) = alpha;

    if opts.display
        fprintf('  %4d   %14.8e   %14.8e   %9.3e  %9.3e  %2d  %6.3f\n', ...
                k, fk, gn, alpha, relres(k+1), bt, model);
    end

    if alpha < opts.alphamin, break; end

    % relative decrease; a long run of tiny decreases with a non-tiny
    % gradient means the search direction is unreliable, not that we have
    % converged -- see the "model" column.
    if (fk - ftry) <= opts.ftol*max(1,abs(fk))
        stag = stag + 1;
        if stag >= opts.stagmax
            exitmsg = sprintf(['stagnated: %d consecutive steps with ' ...
                'relative decrease < %.1e while ||grad||_X = %.3e'], ...
                stag, opts.ftol, gn);
            X = Xtry;
            break
        end
    else
        stag = 0;
    end

    X = Xtry;
end

cost   = cost(1:k+1);
gnorm  = gnorm(1:k+1);
steps  = steps(1:k+1);
relres = relres(1:k+1);

if opts.display
    fprintf('\n  %s\n', exitmsg);
    fprintf('  final cost      : %.8e\n', cost(end));
    fprintf('  final ||grad||_X: %.3e\n', gnorm(end));
    fprintf('  final rel. res  : %.3e\n', relres(end));
end

info = struct('cost',cost,'gnorm',gnorm,'alpha',steps,'relres',relres, ...
              'iter',k,'exitmsg',exitmsg);

% ---------------------------------------------------------------------
% For a correct gradient and small t,
%     (f(X) - f(R_X(-t*g))) / (t*||g||^2)  ->  1.
% A value far from 1 means PROJECTION returned a corrupted direction --
% 
if opts.diagnose
    Gd = hss_euclidean_gradient_mf(X, A);
    gd = projection(Gd);
    q2 = hsstangent_inner(gd, gd);
    fX = f(X);
    fprintf('\n  gradient consistency at the final point:\n');
    fprintf('      t        (f(X)-f(R(-t g)))/(t||g||^2)\n');
    for t = 10.^-(2:6)
        ratio = (fX - f(retraction(X, gd, -t))) / (t*q2);
        fprintf('   %8.1e   %12.6f\n', t, ratio);
    end
    fprintf('   (should approach 1 as t -> 0)\n');
end

% ---------------------------------------------------------------------
% if opts.plots
%     it = 0:k;
%     LW = 2; MS = 5; FS = 12;
% 
%     figure('Color','w','Units','centimeters','Position',[2 2 26 7.5]);
% 
%     subplot(1,3,1);
%     semilogy(it, max(cost,realmin), '-o', 'LineWidth', LW, 'MarkerSize', MS);
%     grid on; box on; set(gca,'FontSize',FS,'LineWidth',1.1);
%     xlabel('iteration'); ylabel('f(X_k)');
%     title('Objective');
% 
%     subplot(1,3,2);
%     semilogy(it, max(gnorm,realmin), '-s', 'LineWidth', LW, 'MarkerSize', MS);
%     grid on; box on; set(gca,'FontSize',FS,'LineWidth',1.1);
%     xlabel('iteration'); ylabel('|| grad f(X_k) ||_X');
%     title('Riemannian gradient norm');
% 
%     subplot(1,3,3);
%     semilogy(it, max(steps,realmin), '-d', 'LineWidth', LW, 'MarkerSize', MS);
%     grid on; box on; set(gca,'FontSize',FS,'LineWidth',1.1);
%     xlabel('iteration'); ylabel('\alpha_k');
%     title('Step size');
% 
%     exportgraphics(gcf,'hss_rgd_history.pdf','ContentType','vector');
% end

end
