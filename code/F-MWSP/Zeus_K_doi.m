function [H,T,C]=Zeus_K_DOI(G)
%generalized set packing toolbox.
%-----------------------------------------
%
%Purpose:  Given G this is the accompanying code for the journal paper.  G describes the entire problem instance and provides links to the necessary sub-routintes which do pricing, and generate DOI.
%
%Input
%
%	G:  Input data structure
%
%		B:  basic info
%			DOI_row_command:  
%			Nd:  number of detections (observations)
%			num_pricing:  number of sub-problems
%			DOI_command:  function name for getting DOI
%			pricing_command:  function name to getting new column				
%			reprice_command:  command for repricing
%			do_sum_list:  set to 1 if we want to add the negative reduced cost of all columns produced in a sub-problem.  USED IN TRACKING BUT NOT LIKELY TO BE USED OTHERWISE
%			pricing_command_heur:  this will be used if it exists AND we have NOT completed an itteration of no returned columns 
%			use_exact:  set to 1 to always use exact.  set to 0 to break when heuristic fails		
%		Q:	ALL user provided computation info for pricing goes here 
%
%
%		opt: options
%
%			max_trip_add:  maximum number of triples to add
%			epsilon:  value for lower bound
%			solve_ilp_prior:  set to 1 if you solve ILP each itteration.  0 otherwise
%			lin_prog_opt:  options for linprog
%			int_lin_prog_opt:  options for intlinprog
%			display_on:   set to 1 if we want to display bound info at each step
%			max_cols_add:  maximum number of columns added
%
%
%			
%
%		T_init:  inital T   OPTIONAL
%
%		C_init
%
%		DOI_options
%
%			num_percentiles to make for each d
%Output 
%
%	H:
%
%		ub:  List of upper bounds  computed across itterations
%		lb:  List of lower bounds computed at each itteration
%		lp:  List of values of the LP relaxation
%		sol:  list of columns produced 
%		aux_info:  auxiliary info produced in each column
%		time_list:  list of times
%			index 1:   master problem 
%			index 2:   pricing problem
%			index 3:   getting the  values of DOI
%			index 4:   producing new triples
%			index 5:   time to solve ILP
%			index 6:   max time of pricing problem
%		pricing_option: list of values across itterationa.
%			value of 0 means heuristic
%			value of 1  mean exact
%			value of 2  means that we have triples 
%		sol_sf:  matrix of association.  has value 1 if two indexs are co-assocated.  on diagonal has 0 if not include in anything
%Important Variables
%		sol_cluster:  vector where index indicates the cluster an observation is associated with.  0 indicates no association 
%
%	T
%		
%		X:    holds along each column the detections found 
%			size:Nd,num_cols
%		Theta:
%			size num_cols,1
%		aux_info
%			size n_info,num_cols
%
%
%	C
%		C3_list(c,k):  k'th value in the triplet c
%			size N_trip by 3
%			
%		C3_mat(ci,d)=1 if d is in triplet ci
%			size N_trip by Nd.  
%
%
%	IMPORTANT
%
%
%	DOI
%		
%		
%		max_val_mat: (d,g) . max amount removing d from any column made up of a subset of g
%
%
%		max_trip_mat: (c,g) max amount removing all but one of the elements could do 
%	
%-----------------------------------------
%structure containing history
H=[];
H.ub=[];
H.lb=[];
H.lp=[];
H.sol=[];
H.aux_info=[];
H.timer=[];
H.pricing_option=[];
%structure containg columns
T=[];%IN THE TEXT THIS IS G,BIG GAMMA
T.X=[];%G_{dg}
T.Theta=[];%GAMMA big 
T.aux_info=[];%auxilliary info for each column
%Structure DOI
DOI=[];
DOI.max_val_mat=[];
DOI.max_trip_mat=[];

%structure containing triplets and doublets
C=[];
C.c3_list=[];
C.c3_mat=[];
if(isfield(G,'C'))
	
	disp('not really setup')
	pause
	C=G.C;
	G.C=[];

end

if(isfield(G,'T'))

	T=G.T;
	G.T=[];
	inds_new=[1:numel(T.Theta)];
	DOI=update_DOI_K(G,C,T,inds_new,DOI);

end

%initialize primal/dual solution
lambda=zeros(G.B.Nd,1);
psi=[];
lp=0;
gamma=[];
%
step=1;

%do heuristic_pricing
pricing_option=0;
%if no pricing command for heuristic is provided turn it off
if(~isfield(G.B,'pricing_command_heur'))
	pricing_option=1;
end
if(~isfield(G.opt,'max_steps'))
	G.opt.max_steps=inf;
end
%START OF THE LOOP FROM 2-8
while(true) 

	%solve the master problem  
	%START OF LINE 3 IN THE CEODE
	t1=tic();
	t_lp=0;
	if(numel(T.Theta)>1.5)%SINCE I MAY INITIALIZE WITH SOME COLUMNS I SOLVE IF THERE ARE COLUMNS PRESENT.  IN THE doc it assumes no columnsare presesent but I may have some.  
		[gamma,lp,lambda,psi,t_lp,LP_SOLVER]=solve_lp_K_DOI(G,T,C,DOI);
		
	
		if(1<0)
			[sol_odd,lp_odd,flag_odd]=linprog(T.Theta,T.X,ones(G.B.Nd,1),[],[],T.Theta*0,(T.Theta*0)+1);
			if(lp-lp_odd>.0001)
				save('oss')
				num_hyp=numel(T.Theta);
				%T.Theta*sol_odd;
				my_sol_candid=[sol_odd;0*LP_SOLVER.slacks(:)];
				gaps=min(LP_SOLVER.B-LP_SOLVER.A*my_sol_candid);
				disp('iossue here')
				lp
				lp_odd
				LP=LP_SOLVER
				[primal_sol_aa,lp_val_aa,flag_aa,~,dual_sol_aa]=cplexlp(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub,[],G.opt.lin_prog_opt);
				[primal_sol_ab,lp_val_ab,flag_ab,~,dual_sol_ab]=linprog(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub,[],G.opt.lin_prog_opt);
				lp_val_aa
				lp_val_ab
				flag_aa
				flag_ab
				save('iss')
				pause
			end
		end
			
	end	
	t1=toc(t1);
	%solve the pricing problem
	t2=tic();
	%lines 4-7 are done internally
	[G,T,residual,inds_new,t6]=do_pricing(G,C,T,lambda,psi,pricing_option);
	t2=toc(t2);
	lb=lp+residual;
	
	%update the DOI the CAPITAL XI_d for all d
	t3=tic();
	DOI=update_DOI_K(G,C,T,inds_new,DOI);
	t3=toc(t3);

	%update list of C
	t4=tic();
	did_add=0;
	if(  residual>-G.opt.epsilon && G.opt.max_trip_add>0.5  &&  (sum((gamma.*(1-gamma)>.01))>0.5) && pricing_option>0.5  )
		[C,did_add]=update_c_list(G,gamma,T,C);
		DOI.max_trip=zeros(size(C.c3_mat,1),1);
		DOI=update_DOI_K(G,C,T,[1:numel(gamma)],DOI);
	end
	t4=toc(t4);

	%Solve the ILP to get rounded solution

	t5=tic();

   	ub=0; 
	cols_keep=[];
	if(G.opt.solve_ilp_prior>0.5)%lines 9-11
		%use the version that reprices ifa vailable
		if(isfield(G.B,'reprice_command'))
			disp('hihih')
			pause
			[ub,cols_keep,T]=solve_ilp_adv(T,G,DOI);
		else
			disp('FIONE')
			%call the version that does not reprice
			[ub,cols_keep]=solve_ilp(T,G);
		end
	end
	t5=toc(t5);
	
	%Store bound info
	H.lb=[H.lb;lb];
	H.lp=[H.lp;lp];
	H.timer=[H.timer;t1,t2,t3,t4,t5,t6,t_lp];
	H.ub=[H.ub;ub];
	H.pricing_option=[H.pricing_option;pricing_option];	
	%store upper bound solution
	[~,ind_min]=min(H.ub);
	if(ind_min==numel(H.ub) && ub<-G.opt.epsilon)
		H.sol=T.X(:,cols_keep);
		H.aux_info=T.aux_info(:,cols_keep);
	end
	%display output
	if(G.opt.display_on)
	
		jy_out_val('step,pricing_option',[step,pricing_option])
		disp('timer info')
		disp('master,pricing,update doi, update c, ilp rounding, max subproblem time')
		jy_out_val('timer',H.timer(end,:))
		disp('note that since I am allowing heuristic lower bound is not valid until termination EVEN when using exat pricing. A subproblem is only solved exactly when the heuristic fails to rpoduce a negative reduced cost solution')
		if(pricing_option<0.5)
		
			jy_out_val('lb,max(lb)',[H.lb(end),max(H.lb)])
		else
			jy_out_val('lb,max(lb(H.pricing_option>0.5))',[H.lb(end),max(H.lb((H.pricing_option>0.5)))])
		end
		jy_out_val('lp',lp)
		jy_out_val('ub,min(ub)',[H.ub(end),min(H.ub)])
		jy_out_val('num_cols,num_rows',[size(T.X,2),size(C.c3_list,1)])
		disp('------')
	end
	%break if ( no triple added AND  no  negative reduced cost found  ) OR tight bound
	cond1=abs(residual)<abs(G.opt.epsilon);
	cond2=did_add<0.5;
	cond3=pricing_option>0.5;
	
	if ( cond1 && cond2 && cond3 )
		disp('breaking') %but im playing time')
	
		if(1==1)
			vals=T.X*gamma;
                        inds_on=find(vals>1.001);
                        my_overage_list=[inds_on,vals(inds_on(:))];
			if(max(vals)>1.001)
				disp('BADS NEWS DOI ACTIVE')
				LP=LP_SOLVER;
				[primal_sol,lp_val,flag,~,dual_sol]=linprog(LP.obj,LP.A,LP.B,[],[],LP.lb,LP.ub,[],G.opt.lin_prog_opt);
				dual_sol=dual_sol.ineqlin;
				lambda=LP.converter.get_lambda*dual_sol;
				gamma=primal_sol(1:numel(T.Theta));
				xi=primal_sol((numel(gamma)+1):end)
				inds_on=find(xi>.0001);
				%LP.rows_use(inds_on,:)
				%jy_out_val('T.Theta(99)',T.Theta(99))
				%this_sol_D=T.X(:,99);

				%this_sol_D(76)=0;
				%my_val_col=G.Q.F.cost.offset+this_sol_D'*G.Q.unary(:)+(0.5*this_sol_D'*G.Q.aff*this_sol_D);				
				%my_val_col
				%DOI
				max(vals)
				save('mr2')
				pause
			end
		end
		break
	end
	if(cond1 && pricing_option<1.5)%update pricing option when usiung heuristic pricing OR IF REGULAR PRICING IS USED BUT NO triples added so far.  IF no violated constraints are found then we will terminate if only heurisitc pricing is used.  
		pricing_option=1;
		if(numel(C.c3_list)>0.5)	
			pricing_option=2;
			disp('triples now here')
			%pause
		else
			disp('pricing is will now use exact solver WHEN heuristic fails')
			%pause
			if(isfield(G.B,'use_exact'))
				if(G.B.use_exact<0.5)
					disp('braking without exact search')
					break
				end
			end
		end
	end
	if(step>=G.opt.max_steps)
		disp('breaking because of too many steps')
		break;
	end
	step=step+1;
end
%solve the ILP rounding problem IF  
if(G.opt.solve_ilp_prior<0.5)%lines 9-11
	ub=[];
	cols_keep=[];
	t5=tic();
	if(isfield(G.B,'reprice_command'))
		[ub,cols_keep,T]=solve_ilp_adv(T,G,DOI);
	else
		[ub,cols_keep]=solve_ilp(T,G);
	end
	t5=toc(t5);
	H.sol=T.X(:,cols_keep);
	H.aux_info=T.aux_info(:,cols_keep);
	H.timer(end,5)=t5;
	H.ub(end)=ub;
	jy_out_val('ub at termination',ub)
end

%COmputing extra formats of solution

v1=[];
v2=[];
n_out=size(H.sol,2);
for(i=1:n_out);
	inds_on=find(H.sol(:,i));
	v1a=repmat(inds_on(:),[numel(inds_on),1]);
	v2a=jy_copy_col(inds_on(:),numel(inds_on));
	v1=[v1;v1a];
	v2=[v2;v2a];
end
v3=(v2*0)+1;
H.sol_sf=sparse(v1,v2,v3,G.B.Nd,G.B.Nd);
%computing second alternative format output
H.sol_cluster=H.sol*([1:size(H.sol,2)]');


disp('done call to Zeus')
