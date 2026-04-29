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

% Tangent fields at this node.
TH.TB12 = alpha * TH0.TB12;
TH.TB21 = alpha * TH0.TB21;

TH.TU = alpha * TH0.TU;
TH.TV = alpha * TH0.TV;
TH.PU = alpha * TH0.PU;
TH.PV = alpha * TH0.PV;

TH.TRl = alpha * TH0.TRl;
TH.PRl = alpha * TH0.PRl;
TH.TRr = alpha * TH0.TRr;
TH.PRr = alpha * TH0.PRr;

TH.TWl = alpha * TH0.TWl;
TH.PWl = alpha * TH0.PWl;
TH.TWr = alpha * TH0.TWr;
TH.PWr = alpha * TH0.PWr;

TH.TD = alpha * TH0.TD;

% Recursive tangent fields.
if ~TH0.leafnode
    TH.TA11 = alpha * TH0.TA11;
    TH.TA22 = alpha * TH0.TA22;
end

end
