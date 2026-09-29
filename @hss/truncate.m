function H = truncate(H, k)
%TRUNCATE Truncate the HSS representation to rank at most K.
%   H = TRUNCATE(H, K) retains up to K singular directions in each
%   nested off-diagonal basis. K must be a finite nonnegative integer.
%   Dense diagonal leaves are unchanged; K does not bound the full matrix rank.
%   Truncation uses SVDs independently of the compression tolerance/options.

validateattributes(k, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'nonnegative'}, 'truncate', 'k');
H = hss_compress(H, double(k), @fixed_rank_svd);

end

function [U, S, V] = fixed_rank_svd(A, k)
[U, S, V] = svd(A, 'econ');
k = min(k, min(size(A)));
U = U(:, 1:k);
S = S(1:k, 1:k);
V = V(:, 1:k);
end
