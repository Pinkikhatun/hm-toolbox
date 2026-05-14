function [TH] = create_randvec(TH, H, r)
%CREATE_RANDVEC Summary of this function goes here
%   Detailed explanation goes here

if H.topnode
    TH = set_hss(TH, H);
end

if ~H.leafnode
    %HSSTANGENT Construct an instance of this class
    %   Detailed explanation goes here
    TH.TB12=randn(r);
    TH.TB21=randn(r);
end
% Factorization of the upper triangular block as U12 * V12'
if ~isempty(H.U)
    %collinear component of U and V
    TH.TU = randn(r);
    TH.TU = TH.TU-TH.TU';
    TH.TV = randn(r);
    TH.TV = TH.TV-TH.TV';
    %%Orthogonal component with respect to U and V
    TH.PU  = randn(size(H, 1),r);
    TH.PU  = TH.PU - TH.U * (TH.U'*TH.PU);
    TH.PV  = randn(size(H, 2),r);
    TH.PV  = TH.PV - TH.V * (TH.V'*TH.PV);
end
% Factorization of the lower triangular block as U21 * V21'
if ~isempty(H.Rl)
    TH.PR = randn(2*r, r);
    TH.PR = TH.PR - [TH.Rl;TH.Rr] * ([TH.Rl;TH.Rr]' * TH.PR);    
    TH.TR = randn(r);
    TH.TR = TH.TR-TH.TR';
    TH.TW = randn(r);
    TH.TW = TH.TW-TH.TW';
    TH.PW = randn(2*r, r);
    TH.PW = TH.PW - [TH.Wl;TH.Wr] * ([TH.Wl;TH.Wr]' * TH.PW);
end
% Dense version of the matrix, if the size is smaller than the
% minimum allowed block size.
if ~isempty(H.D)
    TH.TD = randn(size(H.D));
end

if ~H.leafnode
    TH.TA11 = create_randvec(TH.TA11, H.A11, r);
    TH.TA22 = create_randvec(TH.TA22, H.A22, r);
end


end
