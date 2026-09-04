function TH = create_zerovec(TH, H, r)
%CREATE_ZEROVEC Summary of this function goes here
%   Detailed explanation goes here

if H.topnode
    TH = set_hss(TH, H);
end

if ~H.leafnode
    TH.TB12 = zeros(r);
    TH.TB21 = zeros(r);
end
% Factorization of the upper triangular block as U12 * V12'
if ~isempty(H.U)
    %collinear component of U and V
    TH.TU = zeros(r);
    TH.TV = zeros(r);

    %%Orthogonal component with respect to U and V
    TH.PU  = zeros(size(H.U));
    TH.PV  = zeros(size(H.V));
end

% Factorization of the lower triangular block as U21 * V21'
if ~isempty(H.Rl)
    TH.TR = zeros(r);
    TH.PR = zeros(size([H.Rl; H.Rr]));
    TH.TW = zeros(r);
    TH.PW = zeros(size([H.Wl; H.Wr]));
end
% Dense version of the matrix, if the size is smaller than the
% minimum allowed block size.
if ~isempty(H.D)
    TH.TD = zeros(size(H.D));
end

if ~H.leafnode
    TH.TA11 = create_zerovec(TH.TA11, H.A11, r);
    TH.TA22 = create_zerovec(TH.TA22, H.A22, r);
end


end
