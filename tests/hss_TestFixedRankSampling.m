function hss_TestFixedRankSampling
%HSS_TESTFIXEDRANKSAMPLING Check ranks, sampling widths and approximation.
oldblock = hssoption('block-size');
oldtol = hssoption('threshold');
oldrng = rng;
cleanup = onCleanup(@() restore_options(oldblock, oldtol, oldrng));
hssoption('block-size', 8);
hssoption('threshold', 1e-12);
rng(7);
widths = [];
A = [];
for complex_case = 0:1
    L = randn(32, 3) + complex_case * 1i * randn(32, 3);
    R = randn(24, 3) + complex_case * 1i * randn(24, 3);
    A = L * R';
    for p = [0, 5, 10]
        widths = [];
        args = {3, p};
        if p == 10
            args = {3}; % Default oversampling.
        end
        H = hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), ...
            32, 24, args{:}, 'cluster', [8,16,24,32], [6,12,18,24]);
        assert(isequal(widths, [3+p, 3+p]));
        check_ranks(H, 3);
        assert(norm(full(H) - A, 'fro') < 1e-10 * norm(A, 'fro'));
    end
end
% Rank-deficient input still has exactly the requested basis dimensions.
for scale = [0, 1]
    A = scale * ones(32, 24);
    H = hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), 32, 24, 3);
    check_ranks(H, 3);
    assert(norm(full(H) - A, 'fro') < 1e-10);
end
A = randn(32, 24);
widths = [];
H = hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), 32, 24, 0, 0);
check_ranks(H, 0);
assert(isequal(widths, [0, 0]));
assert(isequal(H.A11.A11.D, A(1:8, 1:6)));
% Fixed rank must not be reduced by the global tolerance.
hssoption('threshold', 1);
H = hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), 32, 24, 3);
check_ranks(H, 3);
hssoption('threshold', 1e-12);
% Existing adaptive path remains usable.
A = randn(32, 2) * randn(2, 32);
H = hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), 32, 32);
assert(norm(full(H) - A, 'fro') < 1e-10 * norm(A, 'fro'));
for args = {{-1}, {1.5}, {Inf}, {NaN}, {[1,2]}, {9}, {3,-1}, {3,1.5}, {3,Inf}}
    failed = false;
    try
        hss('handle', @multiply, @multiply_adjoint, @(i,j) A(i,j), 32, 32, args{1}{:});
    catch
        failed = true;
    end
    assert(failed, 'Invalid rank or oversampling was accepted');
end
fprintf('HSS fixed-rank sampling tests passed.\n');

    function Y = multiply(X)
        widths(end+1) = size(X, 2);
        Y = A * X;
    end
    function Y = multiply_adjoint(X)
        widths(end+1) = size(X, 2);
        Y = A' * X;
    end
end

function check_ranks(H, k)
if H.leafnode
    if ~H.topnode
        assert(size(H.U, 2) == k && size(H.V, 2) == k);
    end
else
    assert(isequal(size(H.B12), [k,k]) && isequal(size(H.B21), [k,k]));
    if ~H.topnode
        assert(isequal(size(H.Rl), [k,k]) && isequal(size(H.Rr), [k,k]));
        assert(isequal(size(H.Wl), [k,k]) && isequal(size(H.Wr), [k,k]));
    end
    check_ranks(H.A11, k);
    check_ranks(H.A22, k);
end
end

function restore_options(block, tol, state)
hssoption('block-size', block);
hssoption('threshold', tol);
rng(state);
end
