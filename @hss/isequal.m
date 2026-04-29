function eq = isequal(H1, H2)
%IS_EQUAL Check if two HSS matrices have the same representation.
%
% EQ = IS_EQUAL(H1, H2) returns true when H1 and H2 have the same HSS tree
% structure and the same data blocks in all nodes.

eq = false;

if ~isa(H1, 'hss') || ~isa(H2, 'hss')
    return;
end

eq = hss_node_equal(H1, H2);

end

function eq = hss_node_equal(H1, H2)
eq = same_size(H1.D, H2.D) && same_values(H1.D, H2.D) && ...
    same_scalar(H1.ml, H2.ml) && same_scalar(H1.nl, H2.nl) && ...
    same_scalar(H1.mr, H2.mr) && same_scalar(H1.nr, H2.nr) && ...
    same_scalar(H1.topnode, H2.topnode) && ...
    same_scalar(H1.leafnode, H2.leafnode) && ...
    same_size(H1.B12, H2.B12) && same_values(H1.B12, H2.B12) && ...
    same_size(H1.B21, H2.B21) && same_values(H1.B21, H2.B21) && ...
    same_size(H1.U, H2.U) && same_values(H1.U, H2.U) && ...
    same_size(H1.V, H2.V) && same_values(H1.V, H2.V) && ...
    same_size(H1.Rl, H2.Rl) && same_values(H1.Rl, H2.Rl) && ...
    same_size(H1.Rr, H2.Rr) && same_values(H1.Rr, H2.Rr) && ...
    same_size(H1.Wl, H2.Wl) && same_values(H1.Wl, H2.Wl) && ...
    same_size(H1.Wr, H2.Wr) && same_values(H1.Wr, H2.Wr);

if ~eq
    return;
end

if H1.leafnode
    eq = isempty(H1.A11) && isempty(H1.A22) && ...
        isempty(H2.A11) && isempty(H2.A22);
else
    eq = isa(H1.A11, 'hss') && isa(H1.A22, 'hss') && ...
        isa(H2.A11, 'hss') && isa(H2.A22, 'hss') && ...
        hss_node_equal(H1.A11, H2.A11) && ...
        hss_node_equal(H1.A22, H2.A22);
end

end

function tf = same_size(A, B)
tf = isequal(size(A), size(B));
end

function tf = same_values(A, B)
if isempty(A) && isempty(B)
    tf = true;
else
    tf = isequal(A, B);
end
end

function tf = same_scalar(a, b)
tf = isequal(a, b);
end

