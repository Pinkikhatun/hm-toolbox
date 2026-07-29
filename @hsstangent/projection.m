function [TH, VH, M, rhs, C] = projection(TH)
%PROJECTION Projection onto the horizontal space.

TH = total_projection(TH);
H = base_point(TH);
r = hssrank(H);
[M, C, rhs] = create_ls_system_rec(TH, r);

%x = M\rhs;
xv = [M'*M, C'; C, zeros(size(C, 1))] \ [M' * rhs; zeros(size(C, 1),1) ];
x = xv(1 : size(M, 2));
%xx = vertical2vec(TH);
%keyboard

VH = create_vertical_from_vec(TH, H, x, r);
tmp = TH - VH;
TH = tmp;
end

function [M, C, rhs] = create_ls_system_rec(TH, r)
    vec = @(X) X(:);

    if TH.leafnode
        M = blkdiag(...
            kron(speye(r), TH.U), ...
            kron(speye(r), TH.V));
        rhs = [ ...
            vec(TH.U * TH.TU) ; ... % + TH.PU - TH.U * (TH.U' * TH.PU)) ; ...
            vec(TH.V * TH.TV) ]; % + TH.PV - TH.V * (TH.V' * TH.PV)) ];
        P = perfect_shuffle(r);        
        C = blkdiag(P, P);
    else
        [M1, C1, rhs1] = create_ls_system_rec(TH.TA11, r);
        [M2, C2, rhs2] = create_ls_system_rec(TH.TA22, r);

        C = blkdiag(C1, C2);

        % We need to add the block corresponding to equations involving C12
        % and C21: these actually refer to the variables Omega and
        % Lambda at the lower level, which are the last 2 r^2 unknowns in
        % the recursive calls. The exact position should be computed
        % looking at the size of M1, M2
        MC = [ ...
            [ sparse(r^2, size(M1, 2)-2*r^2), -kron(TH.B12', speye(r)), sparse(r^2, size(M2,2)), kron(speye(r), TH.B12) ]; ...
            [ sparse(r^2, size(M1, 2)-r^2),   kron(speye(r), TH.B21), sparse(r^2, size(M2,2)-2*r^2), -kron(TH.B21', speye(r)), sparse(r^2, r^2) ]; ...
        ];

        MRW = [...
            [ sparse(r^2, size(M1,2)-2*r^2), -kron(TH.Rl',speye(r)), sparse(r^2, r^2+size(M2,2)) ] ; ...
            [ sparse(r^2, size(M1,2)+size(M2,2)-2*r^2), -kron(TH.Rr',speye(r)), sparse(r^2, r^2) ] ; ...
            [ sparse(r^2, size(M1,2)-r^2), -kron(TH.Wl',speye(r)), sparse(r^2, size(M2,2)) ] ; ...
            [ sparse(r^2, size(M1,2)+size(M2,2)-r^2), -kron(TH.Wr',speye(r)) ] ...
        ];

        M = [ blkdiag(M1, M2); MC ; MRW];
        rhs = [rhs1 ; rhs2 ; vec(TH.TB12) ; vec(TH.TB21) ; ...
            vec(TH.PR(1:r,:) + TH.Rl * TH.TR); 
            vec(TH.PR(r+1:end,:) + TH.Rr * TH.TR); 
            vec(TH.PW(1:r,:) + TH.Wl * TH.TW); 
            vec(TH.PW(r+1:end,:) + TH.Wr * TH.TW); ];

        % Equations for RU * Omega and RV * Lambda
        MRW2 = [ ...
            kron(speye(r), TH.Rl), sparse(r^2, r^2); ...
            kron(speye(r), TH.Rr), sparse(r^2, r^2); ...
            sparse(r^2, r^2), kron(speye(r), TH.Wl); ...
            sparse(r^2, r^2), kron(speye(r), TH.Wr);
        ];

        % Conditions for orthogonality
        P = perfect_shuffle(r);
        C = [ C , sparse(size(C, 1), 2*r^2) ; ...
            [ sparse(r^2, size(M1, 2) - 2*r^2), kron(TH.Rl', TH.Rl'), sparse(r^2, size(M2,2)-r^2), kron(TH.Rr', TH.Rr'), sparse(r^2, 3*r^2) ]; ...
            [ sparse(r^2, size(M1, 2) - r^2), kron(TH.Wl', TH.Wl'), sparse(r^2, size(M2,2)-r^2), kron(TH.Wr', TH.Wr'), sparse(r^2, 2*r^2)]; ...
            [ sparse(2*size(P,1), size(C, 2)) , blkdiag(P,P) ]
        ];
        
        
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
		
		if VH.TA22.leafnode
    		VH.TB12 = VH.B12 * VH.TA22.TV - VH.TA11.TU * VH.B12;
            VH.TB21 = VH.B21 * VH.TA11.TV - VH.TA22.TU * VH.B21;
            VH.PR(:) = [ VH.TA11.TU' * VH.Rl ; VH.TA22.TU' * VH.Rr  ];
            VH.PW(:) = [ VH.TA11.TV' * VH.Wl ; VH.TA22.TV' * VH.Wr  ];
    	else
	    	VH.TB12 = VH.B12 * VH.TA22.TW - VH.TA11.TR * VH.B12;
            VH.TB21 = VH.B21 * VH.TA11.TW - VH.TA22.TR * VH.B21;
            VH.PR(:) = [ VH.TA11.TR' * VH.Rl ; VH.TA22.TR' * VH.Rr  ];
            VH.PW(:) = [ VH.TA11.TW' * VH.Wl ; VH.TA22.TW' * VH.Wr  ];
        end
	end
end

function P = perfect_shuffle(r)
    w = 1 : r^2;
    I = zeros(r, r); I(:) = w;
    I = I';
    I = I(:);
    P = speye(r^2);
    P = P(:, I);
    P = P + speye(r^2);
    J = find(tril(ones(r)) == 1);
    P = P(J, :);
end