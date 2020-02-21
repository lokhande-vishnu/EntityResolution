
function     [gamma,lp_val,lambda,psi,t_lp,LP]=solve_lp(G,T,C,DOI)
%
%Purpose
%
%	solve the master problem given the current primal variables, DOI, and SR-INEQ
%
%Input
%
%
%	G:  Problem instance info
%	T:  Columns
%	C:    SR_Ineq info
%	DOI:  DOI info
%
%Output
%
%	gamma:  primal solution (ignoring xi,xi_3) which are not used
%	lp_val:  objective of primal solution
%	lambda:  dual solution over variables
%	psi: dual solution over triples
%
%------------------------------
%select rows to use 
rows_use=select_rows_use(G,T,C,DOI);
%
LP=form_master_LP_DOI_K(G,T,C,DOI,rows_use);
%LP.obj( (numel(T.Theta)+1):end)=999999;
converter=form_converter_LP_DOI_K(G,T,C,DOI,rows_use);
t_lp=tic();
[primal_sol,lp_val,flag,~,dual_sol]=linprog(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub,[],G.opt.lin_prog_opt);
t_lp=toc(t_lp);
%[j_ilp_val,j_cols_keep]=solve_ilp(T,G);
%j_ilp_val
%lp_val
%save('gogo')
%pause
if(flag<0.5)
	disp('calling cplex')
	[primal_sol,lp_val,flag,~,dual_sol]=cplexlp(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub);%,[],G.opt.lin_prog_opt);

end
if(flag<0.5)
	[primal_sol_cp,lp_val_cp,flag_cp,~,dual_sol_cp]=cplexlp(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub);
	flag_cp
	disp('flag fail')
	flag
	my_gap=T.Theta(:)+LP.A(:,[1:numel(T.Theta)])'*rows_use(:,3);
	jy_out_val('min(my_gap)',min(my_gap));
	save('flagFail')
	pause
end

gamma=primal_sol(1:numel(T.Theta));
lambda=converter.get_lambda*dual_sol.ineqlin;%(1:G.B.Nd);
if(1<0.5)
	
	T.Theta;
	size(lambda)
	size(T.X)
	T.X'*lambda	
	red_cost=T.Theta(:)+(T.X'*lambda);
	red_cost
	if(min(red_cost)<-.0001)
		disp('badNews')
		[j_ilp_val,j_cols_keep]=solve_ilp(T,G);
		save('gogo')
		pause
	end
end
psi=[];
if(numel(C.c3_mat)>0.5)
	psi=converter.get_psi*dual_sol;
end
LP.rows_use=rows_use;
LP.converter=converter;
LP.primal_sol=primal_sol;
LP.slacks=primal_sol(  (1+numel(T.Theta)):end);
