

function  DOI=update_K_DOI(G,C,T,inds_new,DOI)

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
max_val_mat_new=zeros(G.B.Nd,N_new_cols);
num_trip=size(C.c3_list,1);
max_trip_mat_new=zeros(num_trip,N_new_cols);
for(n=1:N_new_cols)%itterate over the columns that are be added 
	
	%get the doi for the specific column
	[max_val,max_trip,DOI]=get_DOI_K_col(G,C,  T.X(  :,inds_new(n)  )  ,  DOI  );
	%update the doi based on this term
	max_val_mat_new(:,n)=max_val;
	if(num_trip>0.5)
		max_trip_mat_new(:,n)=max_trip;
	end
	
	
end

DOI.max_val_mat=[DOI.max_val_mat,max_val_mat_new];
DOI.max_trip_mat=[DOI.max_trip_mat,max_trip_mat_new];

