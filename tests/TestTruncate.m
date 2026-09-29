function TestTruncate(kind)
%TESTTRUNCATE Fixed-rank truncation for HSS and HODLR matrices.
ctor = str2func(kind);
option = str2func([kind 'option']);
oldblock = option('block-size');
oldtol = option('threshold');
cleanup = onCleanup(@() restore_options(option, oldblock, oldtol));
option('block-size', 4);
option('threshold', 1e-14);
oldrng = rng;
rngcleanup = onCleanup(@() rng(oldrng));
rng(42);
for complex_case = 0:1
    A = randn(32, 24) + complex_case * 1i * randn(32, 24);
    H = ctor(A);
    original = full(H);
    for k = [0, 1, 3, 100]
        T = truncate(H, k);
        check_tree(T, H, k);
        assert(isequal(size(T), size(H)));
        if k == 100
            assert(norm(full(T) - original, 'fro') < 1e-11 * norm(original, 'fro'));
        end
    end
    assert(isequal(full(H), original));
    % On a two-leaf tree, each off-diagonal block has its optimal rank-k SVD.
    H = ctor(A, 'cluster', [16, 32], [12, 24]);
    T = truncate(H, 2);
    expected = full(H);
    expected(1:16, 13:24) = best_rank(expected(1:16, 13:24), 2);
    expected(17:32, 1:12) = best_rank(expected(17:32, 1:12), 2);
    assert(norm(full(T) - expected, 'fro') < 1e-11 * norm(expected, 'fro'));
end
for H = {ctor(eye(3)), ctor('zeros', 32, 24), ctor('low-rank', ones(32, 1), ones(24, 1))}
    T = truncate(H{1}, 5);
    check_tree(T, H{1}, 5);
    assert(norm(full(T) - full(H{1}), 'fro') < 1e-11);
end
for k = {-1, 1.5, Inf, NaN, [1, 2], 1i, 'a', []}
    failed = false;
    try
        truncate(ctor(eye(3)), k{1});
    catch
        failed = true;
    end
    assert(failed, 'Invalid rank was accepted');
end
fprintf('%s truncation tests passed.\n', kind);
end

function check_tree(T, H, k)
if isa(T, 'hss')
    leaf = T.leafnode;
    if leaf
        assert(isequal(T.D, H.D));
        if ~T.topnode
            assert(size(T.U, 2) <= k && size(T.V, 2) <= k);
        end
    else
        assert(all(size(T.B12) <= k) && all(size(T.B21) <= k));
        if ~T.topnode
            assert(size(T.Rl, 2) <= k && size(T.Rr, 2) <= k);
            assert(size(T.Wl, 2) <= k && size(T.Wr, 2) <= k);
        end
    end
else
    leaf = isempty(T.A11);
    if leaf
        assert(isequal(T.F, H.F));
    else
        assert(size(T.U12, 2) <= k && size(T.V12, 2) <= k);
        assert(size(T.U21, 2) <= k && size(T.V21, 2) <= k);
    end
end
if ~leaf
    check_tree(T.A11, H.A11, k);
    check_tree(T.A22, H.A22, k);
end
end

function A = best_rank(A, k)
[U, S, V] = svd(A, 'econ');
A = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
end

function restore_options(option, block, tol)
option('block-size', block);
option('threshold', tol);
end
