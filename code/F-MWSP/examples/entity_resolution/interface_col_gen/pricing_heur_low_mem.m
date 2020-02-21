function [heur,tot_resid,new_x_opt,new_x_list,new_aux_list,new_theta_list,new_resid_list,min_col,used_psi,G]=INTERFACE_cell(num_inst,G,C,lambda,psi,this_branch)
%PURPOSE:
%
%	SOlve the pricing problem for cells
%
%Input:
%
%
%	num_inst:  integer stating which is the centroid
%
%	G:  usual
%
%	C:  usual
%
%	psi:  usual
%
%	this_branch:  NEVER WILL BE USED HERE
%
%OUTPUT
%
%	Standard format
%
%		Changes to G. NONE
%-------------------
%hold th
used_psi=[];
new_resid_list=[];
new_theta_list=[];
new_aux_list=[];
tot_resid=[];
new_x_opt=[];
new_x_list=[];

min_col=-1;%note index of min_col
heur=0;
%grab centroid
%get the number of triples
t_ising=tic();

%solve ILP
sol=[];
cost=[];
heur=[];
use_psi=[];
%num_inst
var_keep=find(G.Q.neigh_mat(num_inst,:)>0.5);

anchor=G.Q.anchor_list(num_inst);
if(anchor==0)
	anchor=[];
end
Z=jy_form_ising_bound(G.Q.E,G.Q.pair,G.Q.unary,var_keep,[],0,1,1,G.opt.ILP_opts_cplex,G.opt.int_lin_prog_opt,anchor,-2,G.Q.max_sz_node);
[sol,cost,heur,used_psi]=jy_price_ising_heur(Z,lambda,C.c3_list,psi,G.Q.inst_cost,[],[],G.Q.optimal_ilp_obj(num_inst));

%jy_out_val('[num_inst,toc(t_ising)]',[num_inst,toc(t_ising)])
if(flag<0.5)
	save('flagCase')
	flag
	disp('flag')
	pause
end
sol=sol(1:G.B.Nd);

%add instancing cost 	
%get the solution
if(heur<0)%create one new column if we found a negative reduced cost column
	new_x_opt=sol;
	tot_resid=heur;
	new_x_list=new_x_opt;
        
	%I do not have much to put in the aux term so I just put the centroid
	new_aux_list=[num_inst;0];
	
	%get cost
	new_theta_list=cost;
	new_resid_list=heur;
	min_col=1;
		
end

if(numel(C.c3_list)==0 && sum(lambda)<G.opt.epsilon)
	G.Q.optimal_ilp_obj(num_inst)=heur-G.Q.inst_cost;
end
