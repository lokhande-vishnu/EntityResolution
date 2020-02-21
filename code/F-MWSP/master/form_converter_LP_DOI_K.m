

function converter=form_converter_LP_DOI_K(G,T,C,DOI,rows_use)



	v1=rows_use(:,1);
	v2=[1:numel(v1)];
	v3=(v2*0)+1;
	converter.get_lambda=sparse(v1,v2,v3,G.B.Nd,max(v2));


