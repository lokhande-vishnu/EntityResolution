

function  DOI=update_DOI(G,C,T,inds_new,DOI)

%Purpose
%
%	Update the dual optimal inequalities.  Using the varying DOI 
%	
%
%Input
%
%	G:  usual
%	C:  usual
%	T:  usual
%	inds_new: indexes that correspond to new DOI
%	DOI:  DOI information
%
%Output
%
%	DOI: usual 
%
%
%

N_new_cols=numel(inds_new);
for(n=1:N_new_cols)%itterate over the columns that are be added 
	
	%get the doi for the specific column
	[max_val,max_trip]=get_DOI_col(G,C,  T.X(  :,inds_new(n)  )  ,  T.aux_info(  :,inds_new(n)  )   );
	%update the doi based on this term
	DOI.max_val=(G.opt.epsilon*2)+max([DOI.max_val,max_val],[],2);
	DOI.max_trip=(G.opt.epsilon*2)+max([DOI.max_trip,max_trip],[],2);
	
end



