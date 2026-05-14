function [TH, VH, M, rhs] = projection(TH)
%PROJECTION Projection onto the horizontal space.

TH = total_projection(TH);
H = base_point(TH);
r = hssrank(H);
[M, rhs] = create_ls_system_rec(TH, r);

x = M\rhs;
xx = vertical2vec(TH);
keyboard


VH = create_vertical_from_vec(TH, H, x, r);
tmp = TH - VH;
TH = tmp;
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
        MC = [ ...
            [ sparse(r^2, size(M1, 2)-2*r^2), -kron(TH.B12', speye(r)), sparse(r^2, 2*r^2), kron(speye(r), TH.B12) ]; ...
            [ sparse(r^2, size(M1, 2)-r^2),   kron(speye(r), TH.B21), -kron(TH.B21', speye(r)), sparse(r^2, r^2) ]; ...
        ];

        % Missing: equations local at this level, involving Rl, Rr, Wl, Wr
		%MRW = blkdiag(...
		%  [ sparse(2*r^2, size(M1,2)-2*r^2), blkdiag(-kron(TH.Rl', speye(r)), -kron(TH.Wl', speye(r))) ], ...
		%  [ sparse(2*r^2, size(M2,2)-2*r^2), blkdiag(-kron(TH.Rr', speye(r)), -kron(TH.Wr', speye(r))) ]);
        MRW = [ ...
          [ sparse(r^2, size(M1,2)-2*r^2), -kron(TH.Rl', speye(r)), sparse(r^2, 3*r^2) ]; ...
          [ sparse(r^2, size(M1,2)-r^2), -kron(TH.Wl', speye(r)), sparse(r^2, 2*r^2) ]; ...
		  [ sparse(r^2, size(M1,2)+size(M2,2)-2*r^2), -kron(TH.Rr', speye(r)), sparse(r^2, 1*r^2) ]; ...
          [ sparse(r^2, size(M1,2)+size(M2,2)-r^2), -kron(TH.Wr', speye(r)) ]];

        M = [ blkdiag(M1, M2); MC ; MRW];
        rhs = [rhs1 ; rhs2 ; vec(TH.TB12) ; vec(TH.TB21) ; ...
            vec(TH.PR + [TH.Rl ; TH.Rr ] * TH.TR); vec(TH.PW + [TH.Wl ; TH.Wr ] * TH.TW) ];
        %MRW2 = blkdiag(...
        %    [ kron(speye(r), TH.Rl) ; kron(speye(r), TH.Rr) ], ...
        %    [ kron(speye(r), TH.Wl) ; kron(speye(r), TH.Wr) ]);

        MRW2 = [ blkdiag(kron(speye(r), TH.Rl), kron(speye(r), TH.Wl)); ...
            blkdiag(kron(speye(r), TH.Rr), kron(speye(r), TH.Wr)) ];
        
        
        M = [M, [sparse(size(M1, 1)+size(M2, 1)+2*r^2, size(MRW2, 2)); MRW2]];
        
    end
end
function [VH, l] = create_vertical_from_vec(TH, H, x, r)
	VH = TH;
	VH.TD(:) = 0;
	VH.PV(:) = 0;
	VH.PU(:) = 0;
	VH.PR(:) = 0;
	VH.PW(:) = 0;
	if TH.leafnode
		VH.TU(:) = x(1:r^2);
		VH.TV(:) = x(r^2+ 1:2*r^2);  
		l = 2*r^2;
		VH.TR(:) = 0;
		VH.TW(:) = 0;               
	else    	
		[VH.TA11, l1] = create_vertical_from_vec(TH.TA11, H.A11, x, r);
		[VH.TA22, l2] = create_vertical_from_vec(TH.TA22, H.A22, x(l1+1:end), r);
		
		l = l1 + l2;
		VH.TR(:) = x(l+1:l+r^2); l = l + r^2;
		VH.TW(:) = x(l+1:l+r^2); l = l + r^2;
		VH.TU(:) = 0;
		VH.TV(:) = 0;
		
		if VH.TA11.leafnode
    		VH.TB12 = VH.B12 * VH.TA11.TV - VH.TA22.TU * VH.B12;
    	else
	    	VH.TB12 = VH.B12 * VH.TA11.TW - VH.TA22.TR * VH.B12;
	    end
	    
	    if VH.TA22.leafnode
		    VH.TB21 = VH.B21 * VH.TA22.TV - VH.TA11.TU * VH.B21;
		else
			VH.TB21 = VH.B21 * VH.TA22.TW - VH.TA11.TR * VH.B21;
		end
	end
end
