function     [gamma,lp_val,lambda,psi]=solve_lp(G,T,C,DOI)
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
%create objective
OBJ=[T.Theta;DOI.max_val;DOI.max_trip];
%
A=[T.X];
%Add constraints corresponding to triples
if(numel(C.c3_mat)>0)
	A=[A;(C.c3_mat*T.X)>=1.5];
end
%add in the spot for DOI`
A=[A,-speye(size(A,1))];
%set the upper bound
B=ones(size(A,1),1);
%set the lower/upper bounds
lb=zeros(numel(OBJ),1);
ub=lb+inf;
%grab indexes 
inds_keep=find(OBJ<inf);
%solve master problem 
[primal_sol,lp_val,flag,~,dual_sol]=linprog(OBJ(inds_keep),A(:,inds_keep),B,[],[],lb(inds_keep),ub(inds_keep),[],G.opt.lin_prog_opt);
if(flag<0.5)
	disp('flag fail')
	flag
	save('flagFail')
	%[primal_sol,lp_val,flag,~,dual_sol]=cplexlp(OBJ(inds_keep),A(:,inds_keep),B,[],[],lb(inds_keep),ub(inds_keep),[],G.opt.lin_prog_opt);
	%[primal_sol_1,lp_val,flag,~,dual_sol_1]=linprog(-B,-A',OBJ,[],[],zeros(G.B.Nd,1),inf*ones(G.B.Nd,1),[],G.opt.lin_prog_opt);
	pause
end

%expand the primal solution
primal_sol=sparse(inds_keep,(inds_keep*0)+1,primal_sol,numel(OBJ),1);
%grab primal and dual solution
N_cols=numel(T.Theta);
gamma=primal_sol(1:N_cols);
lambda=dual_sol.ineqlin(1:G.B.Nd);
psi=dual_sol.ineqlin(  (1+G.B.Nd):end);


