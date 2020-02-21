

function [max_val,max_trip,DOI]=entity_doi_K(G,C,this_col,DOI)

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
v3_a=[];
use_original=0;
if(isfield(G.Q,'use_original'));
	use_original=G.Q.use_original;
end
if(use_original>0.5)
	disp('why would I ever want to execute this')
	pause
	v3_a=-sum(G.Q.aff_neg(v1,v1),2);
else
	%disp('not in debujg')
	%pause
	v3_a1=-0.5*sum(G.Q.aff_neg(v1,v1),2);
	v3_a2=-0.5*sum(G.Q.aff(v1,v1),2);
	v3_a=(v3_a1+v3_a2);
	%v3_a(v3_a<0)=0;
end
v3_b=G.Q.unary(v1(:));
v3=.001+v3_a-v3_b;%add in the unary terms
%v3(v3<0)=0;
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



