function TH = mtimes(A, B)
%MTIMES Scalar multiplication for tangent vectors.
%
%   TH = alpha * TH0 or TH = TH0 * alpha
%
% The scalar rescales only tangent components (fields starting with T*),
% while the base-point data is left unchanged.

if isa(A, 'hsstangent') && isscalar(B)
    TH = hsstangent_scalar_mul(B, A);
elseif isa(B, 'hsstangent') && isscalar(A)
    TH = hsstangent_scalar_mul(A, B);
else
    error('hsstangent/mtimes supports only scalar-tangent multiplication');
end

end


function TH = hsstangent_scalar_mul(alpha, TH0)

TH = TH0;

if TH0.leafnode
    TH.TU = alpha * TH0.TU;
    TH.TV = alpha * TH0.TV;
    TH.PU = alpha * TH0.PU;
    TH.PV = alpha * TH0.PV;
    TH.TD = alpha * TH0.TD;
else
    TH.TB12 = alpha * TH0.TB12;
    TH.TB21 = alpha * TH0.TB21;
    TH.TR = alpha * TH0.TR;
    TH.PR = alpha * TH0.PR;
    TH.TW = alpha * TH0.TW;
    TH.PW = alpha * TH0.PW;
    TH.TA11 = alpha * TH0.TA11;
    TH.TA22 = alpha * TH0.TA22;
end

end
