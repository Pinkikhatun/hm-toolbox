function s = hsstangent_inner(TH1, TH2)
%HSSTANGENT_INNER Euclidean inner product in the space of HSS generators.
%
%   s = HSSTANGENT_INNER(TH1, TH2) computes
%
%       s = sum_nodes  <dU1,dU2> + <dV1,dV2> + <dR1,dR2> + <dW1,dW2>
%                    + <dB12_1,dB12_2> + <dB21_1,dB21_2> + <dD1,dD2>
%
%   where each variation is reconstructed from the tangent parameterisation
%
%       dU = U*TU + PU,  dV = V*TV + PV,  dR = R*TR + PR,  dW = W*TW + PW.
%
%   With G = hss_euclidean_gradient(H,A) this gives the directional
%   derivative of f along TH:
%
%       d/dt f( theta + t*TH ) |_{t=0} = hsstangent_inner(G, TH).
%
%   Both arguments must live at the same base point.

s = 0;

if TH1.leafnode
    dU1 = TH1.U * TH1.TU + TH1.PU;
    dU2 = TH2.U * TH2.TU + TH2.PU;
    dV1 = TH1.V * TH1.TV + TH1.PV;
    dV2 = TH2.V * TH2.TV + TH2.PV;

    s = s + dU1(:).' * dU2(:) + dV1(:).' * dV2(:);
    if ~isempty(TH1.TD)
        s = s + TH1.TD(:).' * TH2.TD(:);
    end
else
    R1 = [TH1.Rl; TH1.Rr];   R2 = [TH2.Rl; TH2.Rr];
    W1 = [TH1.Wl; TH1.Wr];   W2 = [TH2.Wl; TH2.Wr];

    if ~isempty(R1)
        dR1 = R1 * TH1.TR + TH1.PR;
        dR2 = R2 * TH2.TR + TH2.PR;
        s = s + dR1(:).' * dR2(:);
    end
    if ~isempty(W1)
        dW1 = W1 * TH1.TW + TH1.PW;
        dW2 = W2 * TH2.TW + TH2.PW;
        s = s + dW1(:).' * dW2(:);
    end

    s = s + TH1.TB12(:).' * TH2.TB12(:);
    s = s + TH1.TB21(:).' * TH2.TB21(:);

    s = s + hsstangent_inner(TH1.TA11, TH2.TA11);
    s = s + hsstangent_inner(TH1.TA22, TH2.TA22);
end

end