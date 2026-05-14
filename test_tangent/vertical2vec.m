function x = vertical2vec(TH)	
	if TH.leafnode
		x = [TH.TU(:); TH.TV(:)];        
	else    	
		x1 = vertical2vec(TH.TA11);
		x2 = vertical2vec(TH.TA22);
		x = [x1; x2; TH.TR(:); TH.TW(:)];
	end
end

