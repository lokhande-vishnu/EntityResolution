function [heur,tot_resid,new_x_opt,new_x_list,new_aux_list,new_theta_list,new_resid_list,ind_opt,used_psi,G]=get_column(m,this_branch,C,lambda,psi,G,pricing_option)
%PURPOSE
%
%	Call the pricing operation for sub-problem m and this branch
%
%INPUT
%
%	m: sub-problem index
%
%	this_branch:  description of the branch
%
%	C:  usual
%
%	lamda:  dual solution over observations
%
%	psi:  dual varibles associated with the triples
%
%	G:  usual	
%	pricing_option:  0 means use heuristic otherwise use exact
%
%OUTPUT
%
%	heur:  heuristic of the pricing problem 
%
%	tot_resid:  this is heur usually.  However in tracking and places where you get many terminations then you get all the lower bounds at each detection
%
%	new_x_opt:  optimal column
%
%	new_x_list:  all columns produced
%
%	new_aux_list:  all aux proeuced
%	
%	new_theta_list:  all cost terms produced
%
%	new_resid_list:  all resid terms produced
%
%	ind_opt:  index of the optimum solution
%
%	used_psi:  list oif indexs that are active in the sub-problem.  This would be the psi terms that the problem is able to incorporate.  
%
%	G:  THis can be updated by the internals

%write command
out_portion='[heur,tot_resid,new_x_opt,new_x_list,new_aux_list,new_theta_list,new_resid_list,ind_opt,used_psi,G]=';
in_portion='(m,G,C,lambda,psi,this_branch);';
fun_portion=G.B.pricing_commands;
if(pricing_option<0.5)
	fun_portion=G.B.pricing_command_heur;
end
my_command=[out_portion,fun_portion,in_portion];


%call evaluation code
eval(my_command);

