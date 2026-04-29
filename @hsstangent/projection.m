function TH = projection(TH)
%PROJECTION Projection onto the horizontal space.

TH = total_projection(TH);

[M, rhs] = create_ls_system_rec(TH, hssrank(base_point(TH)));

end

function [M, rhs] = create_ls_system_rec(TH, r)
    vec = @(X) X(:);

    if TH.leafnode
        M = blkdiag(...
            kron(speye(r), TH.U), ...
            kron(speye(r), TH.V));
        rhs = [ ...
            vec(TH.U * TH.TU) ; ... % + TH.PU - TH.U * (TH.U' * TH.PU)) ; ...
            vec(TH.V * TH.TV) ]; % + TH.PV - TH.V * (TH.V' * TH.PV)) ];
    else
        [M1, rhs1] = create_ls_system_rec(TH.TA11, r);
        [M2, rhs2] = create_ls_system_rec(TH.TA22, r);

        % We need to add the block corresponding to equations involving C12
        % and C21: these actually refer to the variables Omega and
        % Lambda at the lower level, which are the last 2 r^2 unknowns in
        % the recursive calls. The exact position should be computed
        % looking at the size of M1, M2
        MC = [ sparse(r^2, size(M1, 2) - 2*r^2), kron(TH.B12', speye(r)), sparse(r^2, r^2), ...
          sparse(r^2, size(M2, 2) - r^2), kron(speye(r), TH.B12) ; ... 
          sparse(r^2, size(M1, 2) - r^2), kron(TH.B21', speye(r)), kron(speye(r), TH.B21), sparse(r^2, size(M2, 2) - r^2)];

        % Missing: equations local at this level, involving Rl, Rr, Wl, Wr

        M = [ blkdiag(M1, M2); MC ];
        rhs = [rhs1 ; rhs2 ; vec(TH.B12) ; vec(TH.B21) ];
    end
end

