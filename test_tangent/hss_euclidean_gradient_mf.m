function G = hss_euclidean_gradient_mf(H, A)
%HSS_EUCLIDEAN_GRADIENT_MF Euclidean gradient of f(X) = 0.5*||X - A||_F^2.
%
%   G = HSS_EUCLIDEAN_GRADIENT_MF(H, A) returns the gradient of f with
%   respect to every HSS generator of H (leaf bases U,V; transfer matrices
%   Rl,Rr,Wl,Wr; couplings B12,B21; diagonal blocks D), as an HSSTANGENT in
%   the tangent coordinates
%
%       dU = U*TU + PU,  dV = V*TV + PV,
%       dR = R*TR + PR,  dW = W*TW + PW,   R = [Rl;Rr], W = [Wl;Wr],
%       dB12 = TB12,  dB21 = TB21,  dD = TD.

if isa(A, 'hss')
    A = full(A);
end

X = full(H);
if ~isequal(size(X), size(A))
    error('hss_euclidean_gradient_mf:size', ...
        'H and A must have the same size (%dx%d vs %dx%d).', ...
        size(X,1), size(X,2), size(A,1), size(A,2));
end

E = X - A;

G = grad_rec(H, E, [], []);      % [] = the root has no adjoint

end

% -------------------------------------------------------------------------
function G = grad_rec(H, E, PhiU, PhiV)

G = hsstangent();
G.topnode  = H.topnode;
G.leafnode = H.leafnode;

if H.leafnode
    G.U = H.U;   G.V = H.V;   G.D = H.D;

    gU = zero_if_empty(PhiU, size(H.U));   % grad_{U_s} f = Phi_U(s)
    gV = zero_if_empty(PhiV, size(H.V));   % grad_{V_s} f = Phi_V(s)

    G.TU = H.U' * gU;   G.PU = gU - H.U * G.TU;
    G.TV = H.V' * gV;   G.PV = gV - H.V * G.TV;

    G.TD = E;                              % grad_{D_s} f = E(I_s,I_s)

    G.TR = [];  G.PR = [];  G.TW = [];  G.PW = [];
    G.TB12 = []; G.TB21 = [];
else
    m1 = H.ml;  m2 = H.mr;  n1 = H.nl;  n2 = H.nr;
    I1 = 1:m1;          I2 = m1+1 : m1+m2;
    J1 = 1:n1;          J2 = n1+1 : n1+n2;

    E11 = E(I1,J1);  E12 = E(I1,J2);
    E21 = E(I2,J1);  E22 = E(I2,J2);

    G.B12 = H.B12;  G.B21 = H.B21;
    G.Rl  = H.Rl;   G.Rr  = H.Rr;
    G.Wl  = H.Wl;   G.Wr  = H.Wr;
    G.ml  = H.ml;   G.nl  = H.nl;
    G.mr  = H.mr;   G.nr  = H.nr;

    % ---- the four contractions with the children's bases ---------------
    UtE12 = applyUt(H.A11, E12);           % r1 x n2  = U1' * E12
    UtE21 = applyUt(H.A22, E21);           % r2 x n1  = U2' * E21
    E12V2 = applyVt(H.A22, E12.').';       % m1 x r2  = E12 * V2
    E21V1 = applyVt(H.A11, E21.').';       % m2 x r1  = E21 * V1

    % ---- coupling matrices ---------------------------------------------
    G.TB12 = applyUt(H.A11, E12V2);        % U1' * E12 * V2
    G.TB21 = applyUt(H.A22, E21V1);        % U2' * E21 * V1

    % ---- transfer matrices of THIS node --------------------------------
    R = [H.Rl; H.Rr];   W = [H.Wl; H.Wr];

    if isempty(PhiU) || isempty(R)
        G.TR = zeros(size(R,2));  G.PR = zeros(size(R));
    else
        gR   = [ applyUt(H.A11, PhiU(I1,:)) ; applyUt(H.A22, PhiU(I2,:)) ];
        G.TR = R' * gR;           G.PR = gR - R * G.TR;
    end

    if isempty(PhiV) || isempty(W)
        G.TW = zeros(size(W,2));  G.PW = zeros(size(W));
    else
        gW   = [ applyVt(H.A11, PhiV(J1,:)) ; applyVt(H.A22, PhiV(J2,:)) ];
        G.TW = W' * gW;           G.PW = gW - W * G.TW;
    end

    % ---- push the adjoints down to the children ------------------------
    %      local (level-k) term        +  inherited from this node
    PhiU1 = E12V2   * H.B12' + inherit(PhiU, I1, H.Rl);
    PhiU2 = E21V1   * H.B21' + inherit(PhiU, I2, H.Rr);
    PhiV1 = UtE21.' * H.B21  + inherit(PhiV, J1, H.Wl);
    PhiV2 = UtE12.' * H.B12  + inherit(PhiV, J2, H.Wr);

    G.TA11 = grad_rec(H.A11, E11, PhiU1, PhiV1);
    G.TA22 = grad_rec(H.A22, E22, PhiU2, PhiV2);

    G.TU = [];  G.PU = [];  G.TV = [];  G.PV = [];  G.TD = [];
end

end

% -------------------------------------------------------------------------
function Y = inherit(Phi, idx, Rc)
%INHERIT  Phi|_idx * Rc'  -- the term handed down from the parent.
%   Returns the scalar 0 when the parent adjoint is absent (root) or the
%   transfer block is empty, which adds correctly to any matrix and commits
%   to no rank.
if isempty(Phi) || isempty(Rc)
    Y = 0;
else
    Y = Phi(idx,:) * Rc';
end
end

% -------------------------------------------------------------------------
function Y = zero_if_empty(Phi, sz)
if isempty(Phi)
    Y = zeros(sz);
else
    Y = Phi;
end
end

% -------------------------------------------------------------------------
function Y = applyUt(H, X)
%APPLYUT  Y = U^{(k)}' * X, where X has as many rows as the row cluster.
if H.leafnode
    Y = H.U' * X;
else
    Y = H.Rl' * applyUt(H.A11, X(1:H.ml, :)) ...
      + H.Rr' * applyUt(H.A22, X(H.ml+1:end, :));
end
end

% -------------------------------------------------------------------------
function Y = applyVt(H, X)
%APPLYVT  Y = V^{(k)}' * X, where X has as many rows as the column cluster.
if H.leafnode
    Y = H.V' * X;
else
    Y = H.Wl' * applyVt(H.A11, X(1:H.nl, :)) ...
      + H.Wr' * applyVt(H.A22, X(H.nl+1:end, :));
end
end
