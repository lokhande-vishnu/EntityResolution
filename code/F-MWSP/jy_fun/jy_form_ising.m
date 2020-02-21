
function [Z,extra_info]=jy_form_ising(E,pair,unary,var_keep,hard_exclude_cliques,do_repeat_check,has_cplex_opt,has_matlab_opt,cplex_opt,matlab_opt,inds_force_on,inds_integ)
%purpose
%	Solve the ising issue forever
%
%Input
%
%
%	E:  edges as list
%
%	pair:  weights on the edges.  column vector
%
%	unary:  weights on nodes
%
%	var_keep
%		list form: has indexes that we form ising instance over
%	hard_exclude_cliques:  list of cliques that at most one member of can be included
%		this is a matrix with rows descirbing cliques and indexes describing nodes
%
%	has_cplex_opt:  set to 1 if using custom matlab options
%
%	has_mat_opt:  set to 1 if using custom matlab options
%
%	cplex_opt:  cpelx options for ILP.   If nothing provided just display off will be used
%
%	matlab_opt:  matlab options for ILP.  If nothing provided just display off will be used
%	
%	inds_force_on: indexes that must be on
%
%	inds_integ:  defines indexes that are required to be binary
%
%		%leave blank.  use default process of producing them automatically
%		%0:  use fancy LP+rounding proceadure.  NOT done yet.  Ill use 10 second cap on computation
%		%-1:  use none 
%		%-2:   use all single variables
%		%-3:  
%Output
%
%		
%	Z:  complete ising strucutre

%
%		B:  basic info							
%		
%			N_d:  number of detections total NOT JUST in the sub-problem
%			
%			N_e:  number of edges TOTAL not just in sub-problem
%	
%			N_de:  defined to be N_d+N_e
%
%			var_keep:  variables problem is formed over
%	
%			inds_keep:  indexs of the ILP kept in compressed form
%		
%			N_k=numel(inds_keep)
%			
%			E:  same as input
%				
%			unary:  unary terms same as input
%
%			pair: same as input
%	
%			hard_exclude_cliques:  constraints of size Nr, Nd where Nr is the nubme of rwos of hard_exclude_cliques
%				hard_exclud_cliques(r,d)=1 detection d is in row r 
%
%			inds_force_on:  inds to hold on (nodes forced on).  same as input
%
%		M:  matrix conversion info
%			mat_sol_2_full:  takes in the solution in compressed form and converts it to teh complete form
%
%			mat_sol_2_var:  convers the solution in compressed form and converts it to the form of detection
%	
%			mat_dual_2_obj_offset:  takes in dual solution and adds it to the 
%
%		ILP
%
%			OBJ:  Objective
%		
%			OBJ_red:  objective function in compressed form
%	
%			A:  constraint matrix.  only relavant rows included
%
%			A_red:  constraint matrix in compressed form
%			
%			B:  constraint RHS. No need for B_red 
%
%
%			lb:  list of lower bounds 
%
%			lb_red:  lower bound in compressed form		
%
%			ub:  list of upper bounds 
%		
%			ub_red:  upper bound in compressed form
%
%			matlab_integ:  vector of indexes that are required to be integral
%
%			matlab_integ_red:  compressed form of matlab_integ
%
%			cplex_integ:  B,C vector of characters
%
%			cplex_integ_red:  compressed frm of cplex_integ		
%
%		graph
%
%			unary_red
%
%			P_red
%
%
%
%	extra_info:
%	
%		time_gen_pack:  time to do set packing
%-------
%check to see if edges should be corrected.  Re-ordered and 
if(do_repeat_check>0.5)
	%flip incorrectly oriented edges 
	inds_flip=find(E(:,1)>E(:,2));
	E(inds_flip,:)=[E(inds_flip,[2,1])];
	%find non-identical indexes
	inds_keep=find(  E(:,1)~=E(:,2)  );
	E=E(inds_keep,:);
	pair=pair(inds_keep,:);
	
	%retain unique indexes 
	[~,inds_keep]=unique(E,'rows');
	E=E(inds_keep,:);
	pair=pair(inds_keep,:);
end
%

extra_info.time_gen_pack=0;
%identify positive and negative indexes
k1=double(ismember(E(:,1),var_keep));
k2=double(ismember(E(:,2),var_keep));
k3=double(pair<0);
k4=double( pair>0);

%neg indexes are obtained here
neg_inds=find( (k1.*k2.*k3)>0.5  );
neg_inds=neg_inds(:);%make col vector
%pos indexes are obtained 
pos_inds=find( (k1.*k2.*k4)>0.5  );
pos_inds=pos_inds(:);%make col vector 
%store names
N_neg=numel(neg_inds);%number of negative edges relavant
N_pos=numel(pos_inds);%number of positive edes relavant
N_d=size(unary,1);%number of detections
N_e=size(E,1);%number of edges total
N_de=N_d+N_e;%number of nodes and edges 

%construct basic info to Z
Z.B.N_d=N_d;
Z.B.N_e=N_e;
Z.B.N_de=N_de;
Z.B.E=E;
Z.B.var_keep=var_keep;
Z.B.pos_inds=pos_inds;
Z.B.neg_inds=neg_inds;
Z.B.N_neg=numel(neg_inds);
Z.B.N_pos=numel(pos_inds);
Z.B.inds_keep=var_keep(:);
Z.B.inds_force_on=inds_force_on;
if( (numel(neg_inds)+numel(pos_inds) )>0)
	tmp=sort([neg_inds(:);pos_inds(:)]);
	Z.B.inds_keep=[Z.B.inds_keep;(N_d+tmp)];;
end
N_k=numel(Z.B.inds_keep);%numel inds_keep
Z.B.N_k=N_k;%store number of inds_keep
%	form objectiive
%
v1=var_keep(:);
v3=unary(var_keep(:));
%process negative indexes
if(N_neg>0.5)
        v1=[v1;neg_inds(:)+N_d];
        v3=[v3;pair(neg_inds(:))];
end
%process positive indexes
if(N_pos>0.5)
        v1=[v1;pos_inds+N_d];
        v3=[v3;pair(pos_inds)];
end

%write the objective
Z.ILP.OBJ=sparse(v1,(v1*0)+1,v3,N_de,1);%write as a sparse matrix.  Not really sparse but easy ot write
Z.ILP.OBJ_red=Z.ILP.OBJ(Z.B.inds_keep);%hold only the indexes 
A=[];
B=[];
%construct attaractive terms
if(N_neg>0.5)

	v1=jy_copy_col([1:(2*N_neg)  ]',2);
	v3=repmat([-1;1;-1;1],[N_neg,1]);
	
	v2=0*v1;
	v2(1:4:end)=E(neg_inds,1);
	v2(3:4:end)=E(neg_inds,2);
	v2(2:4:end)=N_d+neg_inds;
	v2(4:4:end)=N_d+neg_inds;
	A_neg=sparse(v1,v2,v3,max(v1),N_de);
	A=[A;A_neg];
	%both must be on 
	B=[B;zeros(2*N_neg,1)];

end
%constrcut repulsive edge terms
if(N_pos>0.5)

	v1=jy_copy_col([1:N_pos]',3);
	v2=v1*0;
	v2(1:3:end)=E(pos_inds,1);
	v2(2:3:end)=E(pos_inds,2);	
	v2(3:3:end)=N_d+pos_inds;
	v3=repmat([1;1;-1],[N_pos,1]);
	A_pos=sparse(v1,v2,v3,max(v1),N_de);
        A=[A;A_pos];
	B=[B;ones(N_pos,1)];

end
%construct hard clique constraints 
Z.B.hard_cliques=[];
if(numel(hard_exclude_cliques)>0.5)
	%disp('ttt')
	%jy_out_val('max(sum(hard_cliques,2))',max(sum(hard_exclude_cliques,2)));

	%set indexes to zero for variables not included 
	inds_set_zero=setdiff([1:N_d],var_keep);
	hard_exclude_cliques(:,inds_set_zero)=0;
	%jy_out_val('max(sum(hard_cliques,2))',max(sum(hard_exclude_cliques,2)));

	hard_exclude_cliques(:,inds_force_on)=0;
	%jy_out_val('max(sum(hard_cliques,2))',max(sum(hard_exclude_cliques,2)));

	%identify indexes of rows that should be kept
	rows_keep=find(1.5<sum(hard_exclude_cliques,2));
	%jy_out_val('max(sum(hard_cliques,2))',max(sum(hard_exclude_cliques,2)));
	%add the rows for the hard cliqes
	if(numel(rows_keep)>0)
		A=[A;hard_exclude_cliques(rows_keep,:),zeros(numel(rows_keep),Z.B.N_e)];
		B=[B;(rows_keep*0)+1];
		Z.B.hard_cliques=hard_exclude_cliques(rows_keep,:);
	end
	%disp('fff')
	%pause
end
if(numel(A)==0)%if this is an ising model over one term then we have just that 
	A=zeros(1,Z.B.N_de);
	B=1;
end
%Stores A
Z.ILP.A=A;
Z.ILP.A_red=A(:,Z.B.inds_keep);
%stores B
Z.ILP.B=B;
%store lower bound
Z.ILP.lb=0*Z.ILP.OBJ;
if(numel(inds_force_on)>0.5)%set to one all indexes required to be on which would be anchors
	Z.ILP.lb(inds_force_on)=1;
end

Z.ILP.lb_red=Z.ILP.lb(Z.B.inds_keep);%grab comrpessed terms 
%store upper bound
Z.ILP.ub=1+(0*Z.ILP.OBJ);% grab upper bounds
Z.ILP.ub_red=Z.ILP.ub(Z.B.inds_keep);%grab comrpessed terms 


if(numel(inds_integ)>0)%
	my_command=inds_integ(1);
	if(my_command==0)
		%this will be a fancy LP rounding solution
		disp('not implemented yet')
		jy_mwsp_lp(Z,0);
		pause	
	end
	if(my_command==-3)
		%disp('not implemented yet')
		[inds_integ,extra_info.time_gen_pack]=jy_mwsp_ilp(Z);
		%inds_integ_tmp=jy_get_integ_greedy(Z);
		%jy_out_val('numel(Z.B.var_keep)',numel(Z.B.var_keep))
		%jy_out_val('extra_info.time_gen_pack',extra_info.time_gen_pack)
		%jy_out_val('inds_integ',inds_integ)
		%jy_out_val('inds_integ_tmp',sort(inds_integ_tmp))

		%pause	
	end
	if(my_command==-1)
		%this is nothing is binary required
		inds_integ=[];
	end
	if(my_command==-2)%use all the single variable terms 
		inds_integ=var_keep;
	end	
	
	if(my_command==-5)
		inds_integ=jy_get_integ_greedy(Z);
	end

else	% BIG GREEDY CONSTRUCTION
	disp('bad  but not big deal u need to format input properly')
	pause
	inds_integ=jy_get_integ_greedy(Z);
end

%grab integral indexes for matlab
Z.ILP.matlab_integ=inds_integ;
[~,Z.ILP.matlab_integ_red]=ismember(inds_integ,var_keep);


%set integrality for CPLEX
cplex_integ=67*(ones(1,numel(Z.ILP.OBJ)));
cplex_integ(inds_integ)=66;
Z.ILP.cplex_integ=char(cplex_integ);
Z.ILP.cplex_integ_red=Z.ILP.cplex_integ(Z.B.inds_keep);%matlab_integ_red);
%set cpelx options
if(has_cplex_opt>0.5)%if cplex options provided
	Z.ILP.cplex_opt=cplex_opt;
else
	%set cplex options
	Z.ILP.cplex_opt = cplexoptimset('Display','off');
end
%hold matlab  options
if(has_matlab_opt>0.5)%if matlab options provided
	Z.ILP.matlab_opt=matlab_opt;
else
	%set matlab options
	Z.ILP.matlab_opt=optimset('Display','off');
end

%produce matrix to convert compressed solution to full
v1=Z.B.inds_keep;
v2=[1:Z.B.N_k];
v3=(v2*0)+1;
Z.M.mat_sol_2_full=sparse(v1,v2,v3,N_de,N_k);


%get solution to variales
v1=var_keep;
v2=1:numel(var_keep);
v3=(v2*0)+1;
Z.M.mat_sol_2_var=sparse(v1,v2,v3,N_d,N_k);



%get dual to offset
Z.M.mat_dual_2_obj_offset=sparse([1:numel(var_keep)],var_keep,1+(var_keep*0),N_k,N_d);

%%%%%%%%%
%	graph cut formulation
%
%
%set initial fomula
Z.graph.unary=[];
Z.graph.pairwise_str=[];


if(Z.B.N_pos+Z.B.N_neg+numel(hard_exclude_cliques)>0.5)
	%compute max_val
		
	tmp1=1-1.1*(   sum(unary(Z.B.var_keep)<0));
	
	v1=[Z.B.E(Z.B.neg_inds,1);Z.B.E(Z.B.neg_inds,2)];
	v2=[Z.B.E(Z.B.neg_inds,2);Z.B.E(Z.B.neg_inds,1)];
	v3=pair([Z.B.neg_inds;Z.B.neg_inds]);
	tmp2=max(-sum(sparse(v1,v2,v3,Z.B.N_d,Z.B.N_d),2));
	max_val_ub=1.5*(tmp2+tmp1);
	%tmp2
	%pause;	
	%max_val_ub=1-1.1*(   sum(unary(Z.B.var_keep)<0))-1.1*sum(pair(Z.B.neg_inds));
	
	%v1=pair(:,)
	%save('mr2')
	%max_val_ub
	%pause
	%compute unary terms 
	Z.graph.unary=unary(Z.B.var_keep);
	[~,inds_force_on_red]=ismember(   Z.B.inds_force_on, Z.B.var_keep  );	
	if(numel(inds_force_on_red)>0.5)
		
		Z.graph.unary(inds_force_on_red)=-max_val_ub;
	end
	%compute pairwise terms  for edges 
	pairwise_str=[];
	if(Z.B.N_pos+Z.B.N_neg)
	
	
		pairwise_str=zeros(6,Z.B.N_pos+Z.B.N_neg);%set initial in boykov form
		pairwise_str(6,:)=pair([Z.B.pos_inds(:);Z.B.neg_inds(:)]);%pair terms
		pairwise_str([1,2],:)=[E([Z.B.pos_inds(:);Z.B.neg_inds(:)],:)]';%set indexes
		pairwise_str=sparse(pairwise_str);%make sparse
	end
	
	
	if(numel(hard_exclude_cliques)>0.5)%if hard lciques exist
		for(c=1:size(hard_exclude_cliques,1))%itterate through c
			q=find(hard_exclude_cliques(c,:)>0.5);
			%go through the indexes and add them as pairwise terms
			Nq=numel(q);
			for(i1=1:Nq-1)
				for(i2=i1+1:Nq)
					%add new row for the pairwise matrix	
					pairwise_str=[pairwise_str,[q(i1);q(i2);0;0;0;max_val_ub]];
						
				end
			end
				
			
		
		end

	end
	%project into the reduced formulation
	[~,pairwise_str([1:2],:)]=ismember(   pairwise_str([1:2],:), Z.B.var_keep  );

	%make sparse
	Z.graph.pairwise_str=full( pairwise_str);

end
