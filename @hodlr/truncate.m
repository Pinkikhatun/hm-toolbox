function H = truncate(H, k)
%TRUNCATE Truncate the HODLR representation to rank at most K.
%   H = TRUNCATE(H, K) retains up to K singular directions in each
%   off-diagonal block. K must be a finite nonnegative integer.
%   Dense diagonal leaves are unchanged; K does not bound the full matrix rank.
%   Truncation uses SVDs independently of the compression tolerance/options.

validateattributes(k, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, 'truncate', 'k');
H = truncate_tree(H, double(k));

end

function H = truncate_tree(H, k)
if ~is_leafnode(H)
    H.A11 = truncate_tree(H.A11, k);
    H.A22 = truncate_tree(H.A22, k);
    [H.U12, H.V12] = truncate_factors(H.U12, H.V12, k);
    [H.U21, H.V21] = truncate_factors(H.U21, H.V21, k);
end
end

function [U, V] = truncate_factors(U, V, k)
[QU, RU] = qr(U, 0);
[QV, RV] = qr(V, 0);
[Us, S, Vs] = svd(RU * RV', 'econ');
k = min(k, min(size(S)));
U = QU * Us(:, 1:k) * S(1:k, 1:k);
V = QV * Vs(:, 1:k);
end
