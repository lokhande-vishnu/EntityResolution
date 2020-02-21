function H=Zeus(G)
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
%
%			Nd:  number of detections (observations)
%			num_pricing:  number of sub-problems
%			DOI_command:  function name for getting DOI
%			pricing_command:  function name to getting new column				
%			reprice_command
%			do_sum_list:  set to 1 if we want to add the negative reduced cost of all columns produced in a sub-problem.  USED IN TRACKING BUT NOT LIKELY TO BE USED OTHERWISE
%		
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
%
%
%Important Variables
%
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
%
%	DOI
%		
%		max_val: max_val(d) is the upper bound on dual value d
%			
%
%		max_trip:  max_trip(ci) is the upper bound on dual value for triplet ci
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
%structure containg columns
T=[];
T.X=[];
T.Theta=[];
T.aux_info=[];%auxilliary info for each column

%structure containing triplets and doublets
C=[];
C.c3_list=[];
C.c3_mat=[];

%Structure DOI
DOI=[];
DOI.max_val=[];
DOI.max_trip=[];

%initialize primal/dual solution
lambda=zeros(G.B.Nd,1);
psi=[];
lp=0;
gamma=[];
%
step=1;
while(true) 

	%solve the master problem
	t1=tic();
	if(numel(T.Theta)>1.5)
		[gamma,lp,lambda,psi]=solve_lp(G,T,C,DOI);
	end	
	t1=toc(t1);
	%solve the pricing problem
	%if(numel(C.c3_list)>0.5)
	%	disp('savingtri')
	%	save('savingTir')
	%end
	t2=tic();
	[G,T,residual,inds_new,t6]=do_pricing(G,C,T,lambda,psi);
	t2=toc(t2);

	%update the DOI
	t3=tic();
	DOI=update_DOI(G,C,T,inds_new,DOI);
	t3=toc(t3);
	lb=lp+residual; 

	%update list of C
	t4=tic();
	did_add=0;
	if(  residual>-G.opt.epsilon && G.opt.max_trip_add>0.5  &&  (sum((gamma.*(1-gamma)>.01))>0.5)  )
		[C,did_add]=update_c_list(G,gamma,T,C);
		DOI.max_trip=zeros(size(C.c3_mat,1),1);
		DOI=update_DOI(G,C,T,[1:numel(gamma)],DOI);
	end
	t4=toc(t4);

	%Solve the ILP to get rounded solution

	t5=tic();

   	ub=0; 
	cols_keep=[];
	if(G.opt.solve_ilp_prior>0.5)
		%use the version that reprices ifa vailable
		if(isfield(G.B,'reprice_command'))
			[ub,cols_keep,T]=solve_ilp_adv(T,G,DOI);
		else
			%call the version that does not reprice
			[ub,cols_keep]=solve_ilp(T,G);
		end
	end
	t5=toc(t5);
	
	%Store bound info
	H.lb=[H.lb;lb];
	H.lp=[H.lp;lp];
	H.timer=[H.timer;t1,t2,t3,t4,t5,t6];
	H.ub=[H.ub;ub];
	
	%store upper bound solution
	[~,ind_min]=min(ub);
	if(ind_min==numel(ub) && ub<-G.opt.epsilon)
		H.sol=T.X(:,cols_keep);
		H.aux_info=T.aux_info(:,cols_keep);
	end
	%display output
	if(G.opt.display_on)
	
		jy_out_val('step',step)
		jy_out_val('timer',H.timer(end,:))
		jy_out_val('lb,max(lb)',[H.lb(end),max(H.lb)])
		jy_out_val('lp',lp)
		jy_out_val('ub,min(ub)',[H.ub(end),min(H.ub)])
		jy_out_val('num_cols,num_rows',[size(T.X,2),size(C.c3_list,1)])
		disp('------')
	end
	%break if ( no triple added AND  no  negative reduced cost found  ) OR tight bound
	cond1=abs(residual)<abs(G.opt.epsilon);
	cond2=did_add<0.5;
	cond3=(0==1);%min(H.ub)-max(H.lb)<G.opt.epsilon;

	
	if ( (cond1 && cond2 )|| cond3)
		disp('breaking') %but im playing time')
		break
	end
	step=step+1;
end
%solve the ILP rounding problem IF  
if(G.opt.solve_ilp_prior<0.5)
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
end

