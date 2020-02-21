

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
	for(r=1:G.opt.num_rotations_doi)
		[max_val_rg,max_trip_rg,DOI]=get_DOI_col_rot(G,C,  T.X(  :,inds_new(n)  )   ,DOI,r   );
		%update the doi based on this term
		DOI.max_val_list(:,r)=max([DOI.max_val_list(:,r),max_val_rg],[],2);
		if(numel(C.c3_list)>0.5)
			DOI.max_trip_list(:,r)=max([DOI.max_trip_list(:,r),max_trip_rg],[],2);
		end
	end
end


DOI.max_val=.00001+(sum(DOI.max_val_list,2)/G.opt.num_rotations_doi);


if(numel(C.c3_list)>0.5)

	
	DOI.max_trip=.00001+(sum(DOI.max_trip_list,2)/G.opt.num_rotations_doi);	

end

