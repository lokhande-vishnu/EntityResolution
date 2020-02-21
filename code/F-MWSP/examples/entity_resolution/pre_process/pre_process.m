function G=preprocess_obj_det(F,in_params)
%Purpose
%
%	
%	construct 
%
%
%Input
%
%	F	F is the structure that defines the problem instance
%
%		pairwise: is a |E| by 3 matrix where E is the number of edges (rows in pairwise)
%			pairwise has one row for each pair involved in a cost term
%			pairwise(e,3)=[n1,n2,w_n1_n2].  Here w_n1_n_2 is a real valued weight
%			n_1<n_2 in all cases
%		neigh: is a matrix of E_2 by 2
%			Here E_2 is the number of pairs of edges describing neigborhoods
%			Neigh(e,:)=[n_1,n_2]  indicates that n_1 and n_2 are neighbors. 
			%This means that given that the anchor of a hyptohesis is n_1 then n_2 can be included as a node in the cluster.  
%			Symmetry holds here but does not need to be written
%			n_1<n_2
%		offset:  This is the offset to instance a cluster.  This is theta^0 in the cell paper.
%			default is zero and if not provided I will assume that it is zero.  
%		unary unary which is of size N by 1 where N is the number of nodes
%			unary(n) is the theta_d term in the cell paper (edited)
%			all node indexes are integers in [1:N]
%
%
%	in_params
%		
%
%			
%
%Output
%
%
%	G
%	
%		B
%			Nd.  number of nodes in your graph
%
%
%		Q
%			E:  all edges in pairwise cost terms.  F.pairwise(:,[1,2])
%			pair:  cost terms of pairwise terms.  F.pairwise(:,3);
%			Ne:  number of edges in e.  size(E,1)
%			N:  same as G.B.Nd.  number of nodes in the graph
%			inst_cost: theta^0 term set to zero if not provided
%			aff:  cost terms in pair put into matrix format
%				symmetric but not multiplied by 1/2.  used in calucluaitng cost of new columns during backup AND repricing.  
%			aff_neg:  aff excepte we set all positive tersms to zero.  used in calculating DOI
G=[];%creater structure of problem instance
G.B=[];%create structure of basic info
G.Q=[];%create structure of pricing problem info
F.neigh=double(F.neigh);
G.Q.pair=F.pairwise(:,3);%get cost terms over which pairwise cost terms are described
G.Q.E=F.pairwise(:,[1,2]);%get edges ver wehich pairwise cost terms are described
G.Q.NE=size(G.Q.E,1);%get number of edges

G.B.Nd=max(max(F.pairwise(:,[1,2])));%size(F.pairwise,1);%grab number of detections
G.Q.N=G.B.Nd;%grab number of detections are store in Q

%set inst cost to default if not provided
if(~isfield(F,'inst_cost'))
	F.inst_cost=0;	
end
%store instance cost
G.Q.inst_cost=F.inst_cost;
%construct affinity matrix 
v1=[G.Q.E(:,1);G.Q.E(:,2)];
v2=[G.Q.E(:,2);G.Q.E(:,1)];
v3=[G.Q.pair;G.Q.pair];
G.Q.aff=sparse(v1,v2,v3,G.B.Nd,G.B.Nd);
%set aff_neg as aff except with zero valued entries set to zero.
G.Q.aff_neg=G.Q.aff.*(double(G.Q.aff<0));
%set unary terms to zero if not provided
if(~isfield(F,'unary'))
        F.unary=zeros(G.B.Nd,1);
end
%store unary
G.Q.unary=F.unary;
%-----BIG STEP
%	construct sub-problem list.  
%	We identify every neigborhood that is not a  principal subset of another neighborhood .  

%	step 1.  list all of the neigborhoods 
v1=[F.neigh(:,1);F.neigh(:,2);[1:G.B.Nd]'];
v2=[F.neigh(:,2);F.neigh(:,1);[1:G.B.Nd]'];
v3=(v2*0)+1;
%pretending we ignore unique rows.  
%	neighmat(d_1,d_2)=1 IFF there is is a row in F.neigh of d_1,d_2 OR d_2,d_1  OR d_1=d_2 otherwise zeroi
%a valid column g is one where there exists an anchor d* s.t. G_{dg}=1 implies that d \in Niegh(d_*)
%	d \in Neigh(d_*)
%		IF
%			d=d_*
%			d,d_* is a row in F.neigh
%			d_*,d is a row in F.neigh
%
neigh_mat=sparse(v1,v2,v3);
neigh_mat=sparse(unique(sparse(v1,v2,v3),'rows'));
%step 2 order  the neigborhoods by cardinality.  
[~,my_ord]=sort(sum(neigh_mat,2));

neigh_mat=neigh_mat(my_ord,:);
do_keep=ones(numel(my_ord),1);%hold all non-dominated indexes
num_neigh=numel(my_ord);
for(n=1:num_neigh);	
	%dermine if d is dominated
	row_candid=neigh_mat(n,:);
	row_rem=neigh_mat( (n+1):end,:);
	num_agree=max(row_rem*(row_candid'));
	card_this=sum(row_candid,2);
	if(num_agree>(card_this-0.5))
	%if( max(neigh_mat(d,:)*neigh_mat( (d+1):end,:)) >sum(neigh_mat(d,:))-0.5)
		do_keep(n)=0;
	end
	
end
%neigh(d_1) is dominated by neigh(d_2) if  neigh(d_1) is a PRINCIPAL subst of neigh_d_2
%grab only the neighborhoods 
neigh_mat=neigh_mat(do_keep>0.5,:);
%------------------
%	DONE BIG STEP of computign neighbordhoods
G.B.num_pricing=size(neigh_mat,1);
prob_list=cell(G.B.num_pricing,1);

%store options for ILP sovler
G.opt.int_lin_prog_opt= optimoptions('intlinprog','Display','off');
G.Q.ILP_opts_cplex=cplexoptimset('Display','off');
G.Q.int_lin_prog_opt=optimoptions('intlinprog','Display','off');
G.opt.lin_prog_opt=optimoptions('linprog','Algorithm','interior-point','Display','off');
G.opt.ILP_opts_cplex=cplexoptimset('Display','off');



%initialize time spent solving the cover ILP for ising problem
G.B.tot_time_ilp_cover=0;



for(n=1: G.B.num_pricing  )	
	%grab teh variables that should be kept
	var_keep=find(neigh_mat(n,:)>0.5);
	%construct problem instance using the ising solver
	[prob_list{n},extra_info]=jy_form_ising(G.Q.E,G.Q.pair,G.Q.unary,var_keep,[],0,1,1,G.opt.ILP_opts_cplex,G.opt.int_lin_prog_opt,[],-3);%in_params.opt_get_cover);		
	G.B.tot_time_ilp_cover=G.B.tot_time_ilp_cover+extra_info.time_gen_pack;	
end

%set lower bound on optimal ILP solutions
G.Q.optimal_ilp_obj=-99999999*ones(G.B.num_pricing,1);
G.Q.prob_list=prob_list;
G.B.do_sum_list=1;%means the the lower bound adding the sub-problem costs is used
G.Q.max_doi_bound=9999999999;%set this effectively to infinite
%compute the lwoer bound on what can be acheived in any sub-problem

%provide options to g
%setp max trip add to 1 but can be changed if desired
G.opt.max_trip_add=1;
if(isfield(in_params,'max_trip_add'))
        G.opt.max_trip_add=in_params.max_trip_add;
end

%set to  not solve ilp prior to convervence can be changed if desired
G.opt.solve_ilp_prior=1;
if(isfield(in_params,'solve_ilp_prior'))
	G.opt.solve_ilp_prior=in_params.solve_ilp_prior;
end
%max display on
G.opt.display_on=1;
%set upper bound on the maximum number of columns to add in a given itteariotn.  THis will never be reached
G.opt.max_cols_add=9999;
%provide pricing command
G.B.pricing_commands='pricing_with_trip';
%provide doi command
G.B.DOI_command='obj_det_doi';
%provide epsilon
G.opt.epsilon=.0001;

%set name of command to reprice columns
G.B.reprice_command='reprice_obj'; 

G.Q.force_centroid=0;
