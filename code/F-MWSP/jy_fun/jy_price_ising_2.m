

function [sol,cost,red_cost,used_psi]=jy_price_ising(Z,lambda,C_list,psi,inst_cost,force_on,force_off,lb_obj,group_info)%purpose
%
%	solve pricing problem
%
%
%Input
%
%	Z:  Ising problem instance formulation	
%
%	lambda:  dual variables 
%
%	C_list:  list of triples
%
%	psi
%	
%	inst_cost
%
%	force_on
%
%	force_off
%
%	lb_obj=lower bound on the objective.  
%
%	group_info:  3 indexes
%
%		max_number of groups:  maximum number of groups
%
%		size_groups:  max size of groups
%	
%		skip:  number of indexes to skip between
%
%output
%
%
%	sol:  binary vector having which detections are on and which are off
%	cost:  cost of the solution
%	red_cost:  reduced cost of the solution
%	used_psi:  I flag all triples that are included .  


%if(Z.B.inds_force_on==147)
%	save('kdff')
%end

sol=[];
cost=[];
used_psi=[];

%if(Z.B.inds_force_on==154)
%	disp('moomoo')
%	save('prob')
%end	

if(numel(psi)>0 && 0==1)
	save('dodo')
	disp('I have not debugged or commented that prition yet')
	pause
end	
if(0.5<numel(intersect([Z.B.inds_force_on,force_on],force_off)))%check if the bounds are infeasible 

	disp('BAD')
	pause
	

end

%keep only non-zero valued triples
psi_keep_orig=[];%holds triples that are non-zero valued 
if(numel(psi)>0)
	psi_keep_orig=find(psi>0);
	C_list=C_list(psi>0,:);
	psi=psi(psi>0);
end
use_default=1;%states to use the standard ILP solver (no triples).  THis get 
if(numel(psi)>0.5)
	%NOT DEBUGGED YET
	%disp('not debugged or commented')
	%pause
	
	N_trip=numel(psi);
	v1=jy_copy_col([1:(N_trip*3)]',3);
	v3=repmat([1;1;-1],[N_trip*3,1]);
	v2=v1*0;
	
	v2(1:9:end)=C_list(:,1);
	v2(2:9:end)=C_list(:,2);
	v2(4:9:end)=C_list(:,1);
	v2(5:9:end)=C_list(:,3);
	v2(7:9:end)=C_list(:,2);
	v2(8:9:end)=C_list(:,3);
	v2(3:3:end)=Z.B.N_de+jy_copy_col([1:N_trip]',3);
		
	A_trip_full=sparse(v1,v2,v3,size(C_list,1)*3,Z.B.N_de+N_trip);
	A_trip_compact=[];
	rows_poss=find(1.5<sum(A_trip_full(:,Z.B.var_keep),2));
	psi_keep=[];
	if(numel(rows_poss)>0.5)
		psi_keep=find(-0.5>sum(  A_trip_full(rows_poss, (Z.B.N_de+1):end ) ,1));
		%size(A_trip_full)
		%Z.ILP	
		%save('dongg')
		A_trip_compact=A_trip_full(rows_poss,[Z.B.inds_keep;Z.B.N_de+psi_keep(:)]);
		A=[Z.ILP.A_red,Z.ILP.A_red(:,1)*0;A_trip_compact];
		B=[Z.ILP.B;ones(size(A_trip_compact,1),1)];
		obj=Z.ILP.OBJ_red+Z.M.mat_dual_2_obj_offset*lambda;
		obj=[obj;psi(psi_keep)];
		%is_int=[Z.ilp.is_int,char(67*ones(1,numel(psi)))];
		lb=Z.ILP.lb_red;
        	lb(force_on)=1;
        	ub=Z.ILP.ub_red;
        	ub(force_off)=0;
		lb=[lb;psi_keep*0];
		ub=[ub;1+(psi_keep*0)];
		%[sol_red,red_cost,flag]=cplexmilp(obj,A,B,[],[],[],[],[],Z.ILP.cplex_integ_red,[],Z.ILP.cplex_opt);
		if(~isnan(lb_obj)>0.5)
			A=[A;-obj(:)'];
			B=[B;-lb_obj];
		end
		[~,tmp]=ismember(C_list(:),Z.B.var_keep);
		tmp=tmp(tmp>0.5);
		inds_integ=[Z.ILP.matlab_integ_red(:);tmp(:)];
		inds_integ=unique(inds_integ);
		[sol_red,red_cost,flag]=intlinprog(obj,inds_integ,A,B,[],[],lb,ub,Z.ILP.matlab_opt);

		if(max(sol_red.*(1-sol_red))>.00001)
                	disp('bad')
	
                	disp(max(sol_red.*(1-sol_red)))
                	pause
        	end

		if(flag<0.5)
			disp('flag in triples')
			pause
		end
		used_psi=find(sol_red( (Z.B.N_k+1):end)>0.0001);
		
		if(numel(used_psi)>0)
			used_psi=psi_keep_orig(psi_keep(used_psi));
		end
		sol_red=sol_red(1:Z.B.N_k);
		red_cost=inst_cost+red_cost;
		sol=Z.M.mat_sol_2_full*sol_red;
		cost=inst_cost+(Z.ILP.OBJ(:)'*sol(:));
		use_default=0;
	end
end
if(use_default>0.5)%use_default is triggered if no triples are present or if no triples can be included.  
	%alter bounds 	
	%lb=Z.ILP.lb_red;
	%lb(force_on)=1;
	%ub=Z.ILP.ub_red;
	%ub(force_off)=0;
	obj=Z.ILP.OBJ_red+Z.M.mat_dual_2_obj_offset*lambda;
	%Solve ILP
	sol_red=[];
	red_cost=[];
	flag=1;
	if(numel(obj)==1)
		%[sol_red,red_cost,flag]=cplexmilp(obj,Z.ILP.A_red,Z.ILP.B,[],[],[],[],[],Z.ILP.lb_red,Z.ILP.ub_red,Z.ILP.cplex_integ_red,[],Z.ILP.cplex_opt);
		sol_red=Z.ILP.lb_red+ (double((obj>0))*(Z.ILP.ub_red-Z.ILP.lb_red));
		red_cost=sol_red*obj;
		flag=1;
	else
		A=Z.ILP.A_red;
		B=Z.ILP.B;
		if(~isnan(lb_obj)>0.5)
                        A=[A;-obj(:)'];
			B=[B;-lb_obj];
		end
		%[sol_red,red_cost,flag]=cplexmilp(obj,A,B,[],[],[],[],[],Z.ILP.lb_red,Z.ILP.ub_red,Z.ILP.cplex_integ_red,[],Z.ILP.cplex_opt);
		lb=Z.ILP.lb_red;
		ub=Z.ILP.ub_red;
		A_eq=[];
		B_eq=[];
		if(numel(group_info)>0.5 && numel(Z.ILP.matlab_integ_red)>1.5)
			[A_eq,B_eq]=helper_jy_get_aug_info(Z,group_info);
			n_new_cols=size(Aeq,2)-size(A,2);
			obj=[obj(:);zeros(n_new_cols,1)];
			A=[A,sparse([],[],[],size(A,1),numel(obj_aug) ) ];
			lb=[lb(:);zeros(n_new_cols,1)];
			ub=[ub(:);ones(n_new_cols,1)];
		end		
		[sol_red,red_cost,flag]=intlinprog(obj,Z.ILP.matlab_integ_red,A,B,A_eq,B_eq,lb,ub,Z.ILP.matlab_opt);
	end
	if(max(sol_red.*(1-sol_red))>.00001)
		disp('bad')
		
		disp(max(sol_red.*(1-sol_red)))
		pause
	end	
	if(flag<0.5)
		disp('flag me')
		pause
	end
	%grab solution for returning 
	red_cost=inst_cost+red_cost;
	sol=Z.M.mat_sol_2_full*sol_red;
	cost=inst_cost+(Z.ILP.OBJ'*sol);
end
if(inst_cost<0 && sum(sol)<0.5)%check if the solution is empty and have negative instance cost.  THis means thaty ou have an issue
	disp('you have an issue here.  you cant have negative offset and permissible empty solution')
	pause
end



