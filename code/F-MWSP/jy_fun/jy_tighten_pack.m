

function W=jy_tighen_pack(var_include,X,do_extra_computation)

%Purpose
%
%	generate tightening constraints in a greedy manner
%
%Input
%
%	var_include:  
%
%	X(a,b)=1 if the two variables are not jointly permissable
%	do_extra_computation=1:  if we want to make cliqus of size larger than 2
%Output 
%
%
%	W:  has the cosntraint row
%		W(r,d)=1 if d is part of the constraint r
%
%
%
if(nargin<2.5)
	do_extra_computation=1;	
end
X=double(X);
var_include=var_include(:)';%count variables 
n_inc=numel(var_include);%count the number of varibles 
v1=[];
v2=[];
v3=[];

[~,my_order]=sort(-sum(X(var_include,var_include),2));%cound the number of hard permisavbility cosntriants
my_order=var_include(my_order);
my_order=my_order(:)';%orient my_oder
count=1;%hold the count 
%itterate over pairs 
for(n1=1:(n_inc-1))	
	for(n2=(n1+1):n_inc)

		m1=var_include(n1);
		m2=var_include(n2);
		if(  (X(m1,m2) )<0.5)%if the two are incompatable then we proceed	

			y=X(1,:)*0;%initialize y
			y(m1)=1;%set m1 to be incldued
			y(m2)=1;%set m2 to be included
			if(do_extra_computation>0.5)
				for(n=my_order)%itterate over order from most to least connected 
					
					if ( 0.5>max(X(n,y==1)  ))%it n conflicts with nothing then glag it
						y(n)=1;
					end
				end
			end
			%compute the nodes to include 
			v2b=find(y>0.5);
			v2b=v2b(:);
			%v2b
			if(0==1&& numel(v2b)>2.5)
				disp('dogg')
				v2b
				save('do')
				pause
			end
			%save('dog')
			%add the new 
			v1=[v1;(v2b*0)+count];
			v2=[v2;v2b];
			count=count+1;%increae the count of ht number of rows
	
		end	
	end
end
W=[];
if(numel(v1)>0)
	%Construct w by sparse operation
	W=sparse(v1,v2,(v2*0)+1,max(v1) , size(X,2) );
	
	%get the unique rows
	W=unique(W,'rows');
end

