function Y = hss_perturb(H, TH, t)
%
%   Y = HSS_PERTURB(H, TH, t) applies theta -> theta + t*dtheta to every
%   generator of H, using the tangent parameterisation
%
%       dU = U*TU + PU,  dV = V*TV + PV,  dR = R*TR + PR,  dW = W*TW + PW,
%       dB12 = TB12, dB21 = TB21, dD = TD.

if nargin < 3
    t = 1;
end

Y = H;

if H.leafnode
    Y.U = H.U + t * (H.U * TH.TU + TH.PU);
    Y.V = H.V + t * (H.V * TH.TV + TH.PV);
    Y.D = H.D + t * TH.TD;
else
    Y.A11 = hss_perturb(H.A11, TH.TA11, t);
    Y.A22 = hss_perturb(H.A22, TH.TA22, t);   % NB: H.A22, not H.A11

    R = [H.Rl; H.Rr];
    if ~isempty(R)
        Rnew = R + t * (R * TH.TR + TH.PR);
        m = size(H.Rl, 1);
        Y.Rl = Rnew(1:m, :);
        Y.Rr = Rnew(m+1:end, :);
    end

    W = [H.Wl; H.Wr];
    if ~isempty(W)
        Wnew = W + t * (W * TH.TW + TH.PW);
        m = size(H.Wl, 1);
        Y.Wl = Wnew(1:m, :);
        Y.Wr = Wnew(m+1:end, :);
    end

    Y.B12 = H.B12 + t * TH.TB12;
    Y.B21 = H.B21 + t * TH.TB21;
end

end
