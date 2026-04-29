function TH = set_hss(TH, H)
%SET_HSS Set fields related to the base points of the tangent vector

%Part related HSS matrix
% Factorization of the upper triangular block as U12 * V12'
TH.U = H.U;
TH.V = H.V;

% Factorization of the lower triangular block as U21 * V21'
TH.Rl = H.Rl;
TH.Rr =H.Rr;
TH.Wl = H.Wl;
TH.Wr = H.Wr;

% Size of the matrix
TH.ml =H.ml;
TH.nl =H.nl;
TH.mr =H.mr;
TH.nr =H.nr;

% Dense version of the matrix, if the size is smaller than the
% minimum allowed block size.
TH.D = H.D;

TH.topnode=H.topnode;
TH.leafnode= H.leafnode;

TH.B12 = H.B12;
TH.B21 = H.B21;

if ~H.leafnode
    TH.TA11 = hsstangent;
    TH.TA22 = hsstangent;

    TH.TA11 = set_hss(TH.TA11, H.A11);
    TH.TA22 = set_hss(TH.TA22, H.A22);
end

end

