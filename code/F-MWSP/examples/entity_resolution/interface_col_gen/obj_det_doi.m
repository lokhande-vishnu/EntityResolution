

function [max_val,max_trip]=obj_det_doi(G,C,this_col,this_aux_info)

%Purpose
%
%	Calculate varying DOI
%
%Input
%
%	usual form
%
%OUTput
%
%	usual form
%	

%get indexes for max val (Xi) computation
v1=find(this_col>0.5);
v2=(v1*0)+1;
%add all the pairwise terms
v3=-sum(G.Q.aff_neg(this_col>0.5,this_col>0.5),2);
v3=v3-G.Q.unary(this_col>0.5);%add in the unary terms
%construct matrix.  
max_val=sparse(v1,v2,v3,G.B.Nd,1);
%make that the max value is infinite if thats the center
%initialize triples storage system
num_trip=size(C.c3_list,1);
max_trip=zeros(num_trip,1);
for(c=1:num_trip)       %itterate over triples.  I am lazy so im just using the bound from the smallest two which is sufficeint
	if(1.5<sum(this_col(   C.c3_list(c,:)  )))
		tmp=sort(this_col(C.c3_list(c,:))'.*max_val(C.c3_list(c,:)));
		%grab the smallest two elements
		max_trip(c)=tmp(1)+tmp(2);
                        
	end
end

%if removing the centrodi is ok then add in the negative cost offset
if(G.Q.force_centroid<0.5 && G.Q.inst_cost<0)%
	%NOTE THAT THEIS SHOULD NOT BE PART OF THE COMPUTATION
	disp('dont execute me')
	pause
	max_val=max_val-G.Q.inst_cost;
end


