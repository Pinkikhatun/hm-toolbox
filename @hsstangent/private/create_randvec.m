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
    TH.PU  = randn(size(H, 1)-r,r);
    TH.PV  = randn(size(H, 2)-r,r);
end
% Factorization of the lower triangular block as U21 * V21'
if ~isempty(H.Rl)
    TH.TRl = randn(r);
    TH.TRl = TH.TRl-TH.TRl';
    TH.PRl = randn(H.ml-r, r);
    TH.TRr = randn(r);
    TH.TRr = TH.TRr-TH.TRr';
    TH.PRr = randn(H.mr-r, r);
    TH.TWl = randn(r);
    TH.TWl = TH.TWl-TH.TWl';
    TH.PWl = randn(H.nl-r, r);
    TH.TWr = randn(r);
    TH.TWr = TH.TWr-TH.TWr';
    TH.PWr =randn(H.nr-r, r);
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