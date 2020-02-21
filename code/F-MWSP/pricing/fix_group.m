function new_resid_list=fix_group(new_x_list,new_resid_list,C,psi,used_psi)

%Purpose
%
%	Call fix price on each member of the group.  Specifically we add to new_resid_list the cost terms of triples 
%
%Input
%
%	new_x_list:  list of solutions
%	new_resid_list:  list of residuals assocaited with the new_x_list
%	C:  usual
%	psi:  dual variables of subset row ineuqalities
%	used_psi:  list of psi that have cost already incorporated in the sub-problem cost
%
%Output
%
%	new_resid_list:  the updated resid list 
%

if(numel(psi)>0)%only apply this operation if there is a psi term
	n_sols=size(new_x_list,2);
	for(m=1:n_sols)%itterate over solutions
		x_vec=new_x_list(:,m);%grab the solution
		new_resid_list(m)=fix_price(x_vec,new_resid_list(m),C,psi,used_psi);%adjust the price of the solution

	end
end
