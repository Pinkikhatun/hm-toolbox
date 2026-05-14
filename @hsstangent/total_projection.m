function TH = total_projection(TH)
%TOTAL_PROJECTION Project a tangent vector on the total space.

skew = @(X) .5 *(X - X');

% Projections on the Stiefel manifolds
TH.TU = skew(TH.TU);
TH.TV = skew(TH.TV);
TH.TR = skew(TH.TR);
TH.TW = skew(TH.TW);

if ~TH.leafnode
    TH.TA11 = total_projection(TH.TA11);
    TH.TA22 = total_projection(TH.TA22);
    TH.PR = TH.PR - [TH.Rl;TH.Rr] * ([TH.Rl;TH.Rr]' * TH.PR);
    TH.PW = TH.PW - [TH.Wl;TH.Wr] * ([TH.Wl;TH.Wr]' * TH.PW);
else
    TH.PU  = TH.PU - TH.U * (TH.U'*TH.PU);
    TH.PV  = TH.PV - TH.V * (TH.V'*TH.PV);
end

end

