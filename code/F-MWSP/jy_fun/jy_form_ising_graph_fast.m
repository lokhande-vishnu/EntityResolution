
function [Z,extra_info]=jy_form_ising(E,pair,unary,var_keep,hard_exclude_cliques,conn_mat,do_repeat_check,has_cplex_opt,has_matlab_opt,cplex_opt,matlab_opt,inds_force_on,inds_integ)
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
%Output
%
%		
%	Z:  complete ising strucutre
%
%
%		B:  basic info							
%		
%			N_d:  number of detections total NOT JUST in the sub-problem
%			
%			N_e:  number of edges TOTAL not just in sub-problem
%	
%			N_de:  defined to be N_d+N_e
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
Z.ILP.lb=0*Z.ILP.OBJ;
if(numel(inds_force_on)>0.5)%set to one all indexes required to be on which would be anchors
	Z.ILP.lb(inds_force_on)=1;
end

Z.ILP.lb_red=Z.ILP.lb(Z.B.inds_keep);%grab comrpessed terms 
%store upper bound
Z.ILP.ub=1+(0*Z.ILP.OBJ);% grab upper bounds
Z.ILP.ub_red=Z.ILP.ub(Z.B.inds_keep);%grab comrpessed terms 



%grab integral indexes for matlab


%set integrality for CPLEX
%set cpelx options



%get solution to variales
%get dual to offset
Z.M.mat_dual_2_obj_offset=sparse([1:numel(var_keep)],var_keep,1+(var_keep*0),N_k,N_d);


%%%%%%%%%
%	graph cut formulation
%
%
%set initial fomula
num_keep=size(Z.B.var_keep);
Z.graph.unary=[];
Z.graph.pairwise_str=[];
hard_exclude_pairs=[];
for(n1=1:num_keep-1)
	na=Z.B.var_keep(na);
	nb=Z.B.var_keep(n1+1:num_keep);
	inds_remove=find(conn_mat(n1,nb)==0);
	if(numel(inds_remove>0))
		v2=nb(inds_remove);
		v1=(v2*0)+na;
		hard_exclude_pairs=[hard_exclude_pairs;v1,v2];
	end
end
hard_term=[];
if(numel(hard_exclude_pairs)>0)
	
	n_hard=size(hard_exclude_pairs,1);
	hard_term=[hard_exclude_pairs;zeros(3,n_hard);ones(1,n_hard)];
end	
if(Z.B.N_pos+Z.B.N_neg+numel(hard_exclude_pairs)>0.5)
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
	Z.graph.unary(inds_force_on_red)=-max_val_ub;
	%compute pairwise terms  for edges 
	pairwise_str=[];
	if(Z.B.N_pos+Z.B.N_neg)
	
	
		pairwise_str=zeros(6,Z.B.N_pos+Z.B.N_neg);%set initial in boykov form
		pairwise_str(6,:)=pair([Z.B.pos_inds(:);Z.B.neg_inds(:)]);%pair terms
		pairwise_str([1,2],:)=[E([Z.B.pos_inds(:);Z.B.neg_inds(:)],:)]';%set indexes
		pairwise_str=sparse(pairwise_str);%make sparse
	end
	pairwise_str=[pairwise_str,hard_term];
	if(0==1 && numel(hard_exclude_cliques)>0.5)%if hard lciques exist
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
