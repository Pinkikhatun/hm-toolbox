function M = hssfixedrankfactory(template, rank)
%HSSFIXEDRANKFACTORY Manopt factory for the fixed-rank HSS manifold.
%
%   M = HSSFIXEDRANKFACTORY(H) creates a Manopt manifold whose points are
%   HSS objects with the same tree and (constant) generator rank as H.
%
%   M = HSSFIXEDRANKFACTORY(N, R) first creates a random N-by-N HSS
%   template of rank R.  The current HSS block-size option determines its
%   tree, exactly as for HSSGALLERY('rand', N, R).
%
%   Points are represented by HSS objects and tangent vectors by
%   HSSTANGENT objects.  In particular, an Euclidean gradient supplied to
%   Manopt through problem.egrad must be a gradient with respect to the HSS
%   generators, such as the output of HSS_EUCLIDEAN_GRADIENT_MF.  The
%   factory turns it into a horizontal Riemannian gradient with PROJECTION.
%
%   The factory implements the first-order interface used by Manopt's
%   steepest-descent and conjugate-gradient solvers: metric, projection,
%   retraction, random and zero tangent vectors, linear combinations and a
%   projection-based vector transport.  The geometry and retraction are
%   those implemented by HSSTANGENT_INNER, PROJECTION and RETRACTION.
%
%   Example:
%       H = hssgallery('rand', 16, 2);
%       problem.M = hssfixedrankfactory(H);
%       problem.cost = @(X) .5*norm(full(X) - A, 'fro')^2;
%       problem.egrad = @(X) hss_euclidean_gradient_mf(X, A);
%       X = steepestdescent(problem, H);
%
%   See also HSSTANGENT, HSS_EUCLIDEAN_GRADIENT_MF, STEEPESTDESCENT.

if nargin == 2
    n = template;
    if ~(isnumeric(n) && isscalar(n) && isreal(n) && ...
            isfinite(n) && n == round(n) && n > 0)
        error('hssfixedrankfactory:size', ...
              'N must be a positive integer.');
    end
    if ~(isnumeric(rank) && isscalar(rank) && isreal(rank) && ...
            isfinite(rank) && rank == round(rank) && rank > 0)
        error('hssfixedrankfactory:rank', ...
              'R must be a positive integer.');
    end
    template = hssgallery('rand', n, rank);
elseif nargin ~= 1
    error('hssfixedrankfactory:inputs', ...
          'Use hssfixedrankfactory(H) or hssfixedrankfactory(N, R).');
end

if ~isa(template, 'hss')
    error('hssfixedrankfactory:template', ...
          'The template must be an HSS object.');
end
if template.leafnode
    error('hssfixedrankfactory:singleleaf', ...
          ['The HSS tree must have at least two leaves. Increase N or ' ...
           'decrease hssoption(''block-size'').']);
end

r = size(template.B12, 1);
validate_template(template, r, true);
manifold_dimension = dimension_rec(template, true);
[m, n] = size(template);

M.name = @() sprintf(['Fixed-rank HSS generator quotient manifold ' ...
                      '(%d-by-%d, rank %d)'], m, n, r);
M.dim = @() manifold_dimension;
M.typicaldist = @() sqrt(manifold_dimension);

M.inner = @inner;
M.norm = @norm_tangent;
M.proj = @project;
M.tangent = @project;
M.egrad2rgrad = @project;
M.retr = @retract;
M.rand = @random_point;
M.randvec = @random_tangent;
M.zerovec = @zero_tangent;
M.lincomb = @linear_combination;
M.transp = @transport;

% Tangents already live in the ambient space of HSS generator variations.
M.tangent2ambient_is_identity = true;
M.tangent2ambient = @(~, u) u;

% Small convenience accessors for applications using this factory.
M.rank = @() r;
M.size = @() [m, n];
M.full = @(x) full(x);

    function value = inner(~, u, v)
        require_tangent(u);
        require_tangent(v);
        value = real(hsstangent_inner(u, v));
    end

    function value = norm_tangent(x, u)
        value = sqrt(max(0, inner(x, u, u)));
    end

    function u = project(x, u)
        require_point(x);
        require_tangent(u);
        if ~isequal(base_point(u), x)
            u = rebase_tangent(x, u);
        end
        u = projection(u);
    end

    function y = retract(x, u, t)
        if nargin < 3 || isempty(t)
            t = 1;
        end
        require_point(x);
        require_tangent(u);
        if ~isequal(base_point(u), x)
            error('hssfixedrankfactory:tangentbase', ...
                  'The tangent vector is not based at the supplied point.');
        end
        y = retraction(x, u, t);
    end

    function x = random_point()
        x = randomize_point(template);
    end

    function u = random_tangent(x)
        require_point(x);
        u = random_tangent_rec(x);
        u = project(x, u);
        u_norm = norm_tangent(x, u);
        if ~(isfinite(u_norm) && u_norm > 0)
            error('hssfixedrankfactory:randvec', ...
                  'Failed to generate a nonzero random tangent vector.');
        end
        u = linear_combination(x, 1/u_norm, u);
    end

    function u = zero_tangent(x)
        require_point(x);
        u = zero_tangent_rec(x);
    end

    function u = linear_combination(x, a1, u1, a2, u2)
        require_point(x);
        require_tangent(u1);
        if nargin == 3
            u = tangent_lincomb_rec(x, a1, u1, 0, []);
        elseif nargin == 5
            require_tangent(u2);
            u = tangent_lincomb_rec(x, a1, u1, a2, u2);
        else
            error('hssfixedrankfactory:lincomb', ...
                  'M.lincomb expects either 3 or 5 inputs.');
        end
    end

    function v = transport(~, y, u)
        require_point(y);
        require_tangent(u);
        v = project(y, rebase_tangent(y, u));
    end

    function require_point(x)
        if ~isa(x, 'hss') || ~same_tree(template, x)
            error('hssfixedrankfactory:point', ...
                  'The point must be an HSS object on the factory tree.');
        end
    end

end

% -------------------------------------------------------------------------
function validate_template(H, r, is_root)

if isempty(H.leafnode) || isempty(H.topnode)
    error('hssfixedrankfactory:tree', ...
          'The template has incomplete tree metadata.');
end
if logical(H.topnode) ~= logical(is_root)
    error('hssfixedrankfactory:tree', ...
          'The template has inconsistent topnode flags.');
end

if H.leafnode
    if size(H.U, 2) ~= r || size(H.V, 2) ~= r || ...
            size(H.U, 1) ~= size(H.D, 1) || ...
            size(H.V, 1) ~= size(H.D, 2)
        error('hssfixedrankfactory:nonuniformrank', ...
              ['All leaf bases must have rank %d and agree with their ' ...
               'diagonal block sizes.'], r);
    end
    check_orthonormal(H.U, 'U');
    check_orthonormal(H.V, 'V');
else
    expected = [r, r];
    blocks = {H.B12, H.B21, H.Rl, H.Rr, H.Wl, H.Wr};
    for k = 1:numel(blocks)
        if ~isequal(size(blocks{k}), expected)
            error('hssfixedrankfactory:nonuniformrank', ...
                  ['The current HSSTANGENT implementation requires every ' ...
                   'coupling and transfer block to be %d-by-%d.'], r, r);
        end
    end
    check_orthonormal([H.Rl; H.Rr], 'R');
    check_orthonormal([H.Wl; H.Wr], 'W');
    validate_template(H.A11, r, false);
    validate_template(H.A22, r, false);
end

end

% -------------------------------------------------------------------------
function check_orthonormal(Q, name)

if size(Q, 1) < size(Q, 2)
    error('hssfixedrankfactory:rank', ...
          '%s has more columns than rows.', name);
end
err = norm(Q'*Q - eye(size(Q, 2)), 'fro');
if err > 1e-8*max(1, size(Q, 2))
    error('hssfixedrankfactory:orthonormality', ...
          '%s is not column-orthonormal (error %.3e).', name, err);
end

end

% -------------------------------------------------------------------------
function value = dimension_rec(H, is_root)

if H.leafnode
    ru = size(H.U, 2);
    rv = size(H.V, 2);
    value = numel(H.D) + numel(H.U) - ru^2 + numel(H.V) - rv^2;
else
    R = [H.Rl; H.Rr];
    W = [H.Wl; H.Wr];
    rr = size(R, 2);
    rw = size(W, 2);
    if is_root
        % The root transfer matrices are part of the current HSSTANGENT
        % total-space representation, even though the root has no parent
        % gauge to quotient out.
        dimR = numel(R) - rr*(rr + 1)/2;
        dimW = numel(W) - rw*(rw + 1)/2;
    else
        dimR = numel(R) - rr^2;
        dimW = numel(W) - rw^2;
    end
    value = numel(H.B12) + numel(H.B21) + dimR + dimW + ...
            dimension_rec(H.A11, false) + ...
            dimension_rec(H.A22, false);
end

end

% -------------------------------------------------------------------------
function tf = same_tree(A, B)

tf = logical(A.leafnode) == logical(B.leafnode) && ...
     logical(A.topnode) == logical(B.topnode);
if ~tf
    return;
end

if A.leafnode
    tf = isequal(size(A.D), size(B.D)) && ...
         isequal(size(A.U), size(B.U)) && ...
         isequal(size(A.V), size(B.V));
else
    tf = isequal([A.ml, A.mr, A.nl, A.nr], ...
                 [B.ml, B.mr, B.nl, B.nr]) && ...
         isequal(size(A.B12), size(B.B12)) && ...
         isequal(size(A.B21), size(B.B21)) && ...
         isequal(size(A.Rl), size(B.Rl)) && ...
         isequal(size(A.Rr), size(B.Rr)) && ...
         isequal(size(A.Wl), size(B.Wl)) && ...
         isequal(size(A.Wr), size(B.Wr)) && ...
         same_tree(A.A11, B.A11) && same_tree(A.A22, B.A22);
end

end

% -------------------------------------------------------------------------
function H = randomize_point(H)

if H.leafnode
    H.U = random_stiefel(size(H.U, 1), size(H.U, 2));
    H.V = random_stiefel(size(H.V, 1), size(H.V, 2));
    H.D = randn(size(H.D));
else
    H.B12 = randn(size(H.B12));
    H.B21 = randn(size(H.B21));

    R = random_stiefel(size(H.Rl, 1) + size(H.Rr, 1), size(H.Rl, 2));
    split = size(H.Rl, 1);
    H.Rl = R(1:split, :);
    H.Rr = R(split+1:end, :);

    W = random_stiefel(size(H.Wl, 1) + size(H.Wr, 1), size(H.Wl, 2));
    split = size(H.Wl, 1);
    H.Wl = W(1:split, :);
    H.Wr = W(split+1:end, :);

    H.A11 = randomize_point(H.A11);
    H.A22 = randomize_point(H.A22);
end

end

% -------------------------------------------------------------------------
function Q = random_stiefel(m, r)

[Q, R] = qr(randn(m, r), 0);
signs = sign(diag(R));
signs(signs == 0) = 1;
Q = Q*diag(signs);

end

% -------------------------------------------------------------------------
function u = random_tangent_rec(H)

u = zero_tangent_rec(H);
if H.leafnode
    [u.TU, u.PU] = random_stiefel_tangent(H.U);
    [u.TV, u.PV] = random_stiefel_tangent(H.V);
    u.TD = randn(size(H.D));
else
    u.TB12 = randn(size(H.B12));
    u.TB21 = randn(size(H.B21));
    [u.TR, u.PR] = random_stiefel_tangent([H.Rl; H.Rr]);
    [u.TW, u.PW] = random_stiefel_tangent([H.Wl; H.Wr]);
    u.TA11 = random_tangent_rec(H.A11);
    u.TA22 = random_tangent_rec(H.A22);
end

end

% -------------------------------------------------------------------------
function [TQ, PQ] = random_stiefel_tangent(Q)

Z = randn(size(Q));
dQ = Z - Q*((Q'*Z + Z'*Q)/2);
TQ = Q'*dQ;
PQ = dQ - Q*TQ;

end

% -------------------------------------------------------------------------
function u = zero_tangent_rec(H)

u = hsstangent();
u.B12 = H.B12;
u.B21 = H.B21;
u.U = H.U;
u.V = H.V;
u.Rl = H.Rl;
u.Rr = H.Rr;
u.Wl = H.Wl;
u.Wr = H.Wr;
u.ml = H.ml;
u.nl = H.nl;
u.mr = H.mr;
u.nr = H.nr;
u.D = H.D;
u.topnode = H.topnode;
u.leafnode = H.leafnode;

if H.leafnode
    u.TU = zeros(size(H.U, 2));
    u.TV = zeros(size(H.V, 2));
    u.PU = zeros(size(H.U));
    u.PV = zeros(size(H.V));
    u.TD = zeros(size(H.D));
else
    u.TB12 = zeros(size(H.B12));
    u.TB21 = zeros(size(H.B21));
    u.TR = zeros(size(H.Rl, 2));
    u.PR = zeros(size([H.Rl; H.Rr]));
    u.TW = zeros(size(H.Wl, 2));
    u.PW = zeros(size([H.Wl; H.Wr]));
    u.TA11 = zero_tangent_rec(H.A11);
    u.TA22 = zero_tangent_rec(H.A22);
end

end

% -------------------------------------------------------------------------
function u = tangent_lincomb_rec(H, a1, u1, a2, u2)

u = zero_tangent_rec(H);
if H.leafnode
    u.TU = combine(a1, u1.TU, a2, u2, 'TU');
    u.TV = combine(a1, u1.TV, a2, u2, 'TV');
    u.PU = combine(a1, u1.PU, a2, u2, 'PU');
    u.PV = combine(a1, u1.PV, a2, u2, 'PV');
    u.TD = combine(a1, u1.TD, a2, u2, 'TD');
else
    u.TB12 = combine(a1, u1.TB12, a2, u2, 'TB12');
    u.TB21 = combine(a1, u1.TB21, a2, u2, 'TB21');
    u.TR = combine(a1, u1.TR, a2, u2, 'TR');
    u.PR = combine(a1, u1.PR, a2, u2, 'PR');
    u.TW = combine(a1, u1.TW, a2, u2, 'TW');
    u.PW = combine(a1, u1.PW, a2, u2, 'PW');
    if isempty(u2)
        u.TA11 = tangent_lincomb_rec(H.A11, a1, u1.TA11, 0, []);
        u.TA22 = tangent_lincomb_rec(H.A22, a1, u1.TA22, 0, []);
    else
        u.TA11 = tangent_lincomb_rec(H.A11, a1, u1.TA11, a2, u2.TA11);
        u.TA22 = tangent_lincomb_rec(H.A22, a1, u1.TA22, a2, u2.TA22);
    end
end

end

% -------------------------------------------------------------------------
function value = combine(a1, value1, a2, u2, field)

if isempty(u2)
    value = a1*value1;
else
    value = a1*value1 + a2*u2.(field);
end

end

% -------------------------------------------------------------------------
function v = rebase_tangent(H, u)

v = zero_tangent_rec(H);
if H.leafnode
    dU = u.U*u.TU + u.PU;
    dV = u.V*u.TV + u.PV;
    v.TU = H.U'*dU;
    v.PU = dU - H.U*v.TU;
    v.TV = H.V'*dV;
    v.PV = dV - H.V*v.TV;
    v.TD = u.TD;
else
    oldR = [u.Rl; u.Rr];
    newR = [H.Rl; H.Rr];
    dR = oldR*u.TR + u.PR;
    v.TR = newR'*dR;
    v.PR = dR - newR*v.TR;

    oldW = [u.Wl; u.Wr];
    newW = [H.Wl; H.Wr];
    dW = oldW*u.TW + u.PW;
    v.TW = newW'*dW;
    v.PW = dW - newW*v.TW;

    v.TB12 = u.TB12;
    v.TB21 = u.TB21;
    v.TA11 = rebase_tangent(H.A11, u.TA11);
    v.TA22 = rebase_tangent(H.A22, u.TA22);
end

end

% -------------------------------------------------------------------------
function require_tangent(u)

if ~isa(u, 'hsstangent')
    error('hssfixedrankfactory:tangent', ...
          'Tangent vectors must be HSSTANGENT objects.');
end

end
