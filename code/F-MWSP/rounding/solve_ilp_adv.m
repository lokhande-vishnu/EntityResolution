function     [ilp_val,cols_keep,T]=solve_ilp(T,G,DOI)
%Purpose
%
%	Produce a rounded solution and new columns. THe proceadure is that we itterate between solving hte ILP and addding new columns until convergence
%
%Input
%
%	T: Holds current columns
%	G: holds problem structure
%
%Output
%
%	ilp_val:  Objective of integer solution
%	cols_keep:  grabs the columns selected in the ILP solution
%	T:  updated with the new cols
%Formulate ILP
did_term=0;
ilp_val=[];
cols_keep=[];
while(true)
	%se up ILP
	OBJ=[T.Theta;DOI.max_val(:)];%- is because matlab minimizes not maximizes the value of LPs
	N_hyp=numel(T.Theta);
	A=[T.X,-speye(G.B.Nd)];
	B=ones(size(A,1),1);
	lb=0*OBJ;
	ub=lb+inf;
	%Solve ILP
	%ones(1,N_hyp)
	bin_vec=67+(OBJ*0);
	bin_vec(1:N_hyp)=66;
	bin_vec=char(bin_vec(:)');
	[primal_sol,ilp_val]=cplexmilp(OBJ,A,B,[],[],[],[],[],lb,ub,bin_vec,[],G.Q.ILP_opts_cplex);
	%[primal_sol,ilp_val]=intlinprog(OBJ,[1:N_hyp],A,B,[],[],lb,ub,G.opt.int_lin_prog_opt );
	
	%grab solution
	cols_keep=find(primal_sol(1:N_hyp)>0.5);
	Xi=primal_sol((1+N_hyp):end);
	
	%break if no Xi available
	if(sum(Xi)<0.5)
		did_term=1;
		break;
	end
	
	%find rows to set to zero
	rows_kill=find(Xi>0.5);
	cols_keep=cols_keep(:)';%cols_keep(1,:);
	for(k=1:numel(cols_keep))
		%grab column and set its values on killed indexes to zero and call repreice
		candid_col=T.X(:,cols_keep(k));
		if(sum(candid_col(rows_kill))>0.0001)%operate if there are indexes to kill
			candid_col(rows_kill)=0;
			T2=get_reprice(G,candid_col);	%reprice it
		
			if(numel(T2.Theta)>0)%if there are elements that have cost that are non-negative
				inds_add=ismember(T2.X',T.X','rows');
				inds_add=find(inds_add<0.5);
				inds_add=inds_add(:)';
				if(numel(inds_add)>0.5)%add those elements
					T.X=[T.X,T2.X(:,inds_add)];
					T.Theta=[T.Theta;T2.Theta(inds_add)];
					T.aux_info=[T.aux_info,T2.aux_info(:,inds_add)];
				end
			end
		end
	end	
end

