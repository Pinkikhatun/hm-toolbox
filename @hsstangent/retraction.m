function Y = retraction(H, TH, t)
%RETRACTION Retraction on the HSS quotient manifold.
%
%   Y = RETRACTION(H, TH, t) moves the base point H along the tangent vector
%   TH by the step t, re-orthonormalising every factor with UF (the polar /
%   SVD factor), which is what keeps the iterate on the manifold.
%
%   TWO CORRECTIONS with respect to the previous version:
%     (1) "if nargin < 2" was testing the wrong argument, so calling
%         RETRACTION(H,TH) with two inputs left t undefined. It must be
%         "nargin < 3".
%     (2) the recursive call for the second child read
%             Y.A22 = retraction(H.A11, TH.TA22, t);
%         and has been changed to H.A22. With H.A11 the right subtree was
%         retracted from the WRONG base point, which corrupts every iterate
%         of a gradient-descent loop without raising an error.

if nargin < 3 || isempty(t)
    t = 1;
end

Y = H;

if H.leafnode
    %% Retract U
    dU  = H.U * TH.TU + TH.PU;
    Y.U = uf(H.U + t*dU);

    %% Retract V
    dV  = H.V * TH.TV + TH.PV;
    Y.V = uf(H.V + t*dV);

    %% Dense block
    Y.D = H.D + t*TH.TD;
else
    %% Recursive HSS blocks
    Y.A11 = retraction(H.A11, TH.TA11, t);
    Y.A22 = retraction(H.A22, TH.TA22, t);      % was H.A11

    %% Retract transfer matrix R
    R    = [H.Rl; H.Rr];
    dR   = R*TH.TR + TH.PR;
    Rnew = uf(R + t*dR);
    m     = size(H.Rl,1);
    Y.Rl  = Rnew(1:m,:);
    Y.Rr  = Rnew(m+1:end,:);

    %% Retract transfer matrix W
    W    = [H.Wl; H.Wr];
    dW   = W*TH.TW + TH.PW;
    Wnew = uf(W + t*dW);
    m     = size(H.Wl,1);
    Y.Wl  = Wnew(1:m,:);
    Y.Wr  = Wnew(m+1:end,:);

    %% Core blocks
    Y.B12 = H.B12 + t*TH.TB12;
    Y.B21 = H.B21 + t*TH.TB21;
end

end

function A = uf(A)
    [L, ~, R] = svd(A, 0);
    A = L*R';
end
