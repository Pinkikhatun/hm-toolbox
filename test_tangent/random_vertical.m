function [VH, l] = random_vertical(TH, r)	
	VH = TH;
	VH.TD(:) = 0;
	VH.PV(:) = 0;
	VH.PU(:) = 0;
	VH.PR(:) = 0;
	VH.PW(:) = 0;
	if TH.leafnode
		VH.TU(:) = randn(r^2, 1);
		VH.TU(:) = skew(VH.TU);
		VH.TV(:) = randn(r^2, 1);  
		VH.TV(:) = skew(VH.TV);
		l = 2*r^2;
		VH.TR(:) = 0;
		VH.TW(:) = 0;               
	else    	
		[VH.TA11, l1] = random_vertical(TH.TA11, r);
		[VH.TA22, l2] = random_vertical(TH.TA22, r);
		
		l = l1 + l2;
		VH.TR(:) = randn(r^2, 1); l = l + r^2;
		VH.TR(:) = skew(VH.TR);		
		VH.TW(:) = randn(r^2, 1); l = l + r^2;
		VH.TW(:) = skew(VH.TW);
		VH.TU(:) = 0;
		VH.TV(:) = 0;
		
		if VH.TA11.leafnode
    		VH.TB12 = VH.B12 * VH.TA22.TV - VH.TA11.TU * VH.B12;
    	else
	    	VH.TB12 = VH.B12 * VH.TA22.TW - VH.TA11.TR * VH.B12;
	    end
	    
	    if VH.TA22.leafnode
		    VH.TB21 = VH.B21 * VH.TA11.TV - VH.TA22.TU * VH.B21;
		else
			VH.TB21 = VH.B21 * VH.TA11.TW - VH.TA22.TR * VH.B21;
		end
	end
end

function A = skew(B)
 	A = 0.5*(B-B');
end
