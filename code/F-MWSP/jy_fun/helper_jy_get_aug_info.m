function [A_eq,B_eq]=helper_jy_get_aug_info(Z,group_info)



Aeq=[];

my_groups=[];

for(j=1:group_info(3))

	start_index=1+(group_info(2)*(j-1));
	stop_index=min(start_index-1+group_info(1),numel(Z.ILP.matlab_integ_red));
	this_group={[Z.ILP.matlab_integ_red(start_index:stop_index)]};
	my_groups=[my_groups;this_group];
	if(stop_index==numel(Z.ILP.matlab_integ_red))
	
		break
	end	
	
end
B_eq=[];
for(j=1:numel(my_groups))


	
	n_var_j=numel(my_groups{j}(:)));
	A_eq_new=jy_make_bl	
	v3a_1=repmat([-1;1]);
	

	v3a=[v3a_1;v3a_2];
	if(j>1.5)
		v1a=v1a+max(v1);
	end
	v1=[v1;v1a];
	v2=[v2;v2a];
	v3=[v3;v3a];
	n_new_cols	
	
	A_eq=[A_eq,sparse([],[],[],size(A_eq,1),n_new_cols);A_eq_new]	
	B_eq=[B_eq;zeros(n_var_j,1);1];

end







