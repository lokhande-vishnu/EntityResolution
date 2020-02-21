

function [sol,cost,red_cost,used_psi]=jy_price_ising_heur(Z,lambda,C_list,psi,inst_cost,force_on,force_off,lb_obj)
%purpose
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
%output
%
%
%	sol:  binary vector having which detections are on and which are off
%	cost:  cost of the solution
%	red_cost:  reduced cost of the solution
%	used_psi:  I flag all triples that are included .  


sol=[];
cost=[];
used_psi=[];

%keep only non-zero valued triples
obj=Z.ILP.OBJ_red+Z.M.mat_dual_2_obj_offset*lambda;
%Solve ILP
sol_red=[];
red_cost=[];
flag=1;
if(numel(obj)==1 ||  numel(Z.graph.unary)<0.5)
	sol_red=Z.ILP.lb_red+ (double((obj>0))*(Z.ILP.ub_red-Z.ILP.lb_red));
	red_cost=sol_red*obj;
	
	flag=1;

	red_cost=inst_cost+red_cost;
	sol=Z.M.mat_sol_2_full*sol_red;
	cost=inst_cost+(Z.ILP.OBJ'*sol);

end
if(numel(obj)>1)

	A=Z.ILP.A_red;
	B=Z.ILP.B;
	if(~isnan(lb_obj)>0.5)
		A=[A;-obj(:)'];
		B=[B;-lb_obj];
	end
	lb=Z.ILP.lb_red;
	ub=Z.ILP.ub_red;
	use_heur=1;
	unary_w_dual=[lambda(Z.B.var_keep)*0,Z.graph.unary+lambda(Z.B.var_keep)];
	%unary_w_dual(Z.B.,1);
	gc_sol=int32(Z.B.var_keep*0);
	try
		if(0==1)
			%if (Z.B.inds_force_on==107)
			%	save('look107')	
				%end
			%disp('probing')
			max_abs_pair=max(max(abs(Z.graph.pairwise_str(3:6,:))));
			max_abs_unary=max(abs(unary_w_dual(:)));
			if(max_abs_pair>100 || max_abs_unary>100  )
				disp('myGoose')
				save('myGoose')
				pause
			end
			save('mym')
			gc_sol = QPBO_wrapper_mex(unary_w_dual', Z.graph.pairwise_str, int32(Z.B.var_keep*0), 'p');
		end
	catch
		disp('error in probe')
		save('probeError')
	end
	try
		%disp('itteating ')
		gc_sol = QPBO_wrapper_mex(unary_w_dual', Z.graph.pairwise_str, gc_sol, 'i');
	catch
		disp('error in improve')
		save('errorImprove')
		pause
	end
	gc_sol=double(gc_sol);
	lb(  gc_sol==1 )=1;
	ub(gc_sol==0)=0;
	ub(gc_sol<0.5)=0;
	gc_sol(gc_sol<0.5)=0;	
	%cost_1=sum(lambda(Z.B.var_keep))
	%xi1=gc_sol(Z.graph.pairwise_str(1,:));
	%xi2=gc_sol(Z.graph.pairwise_str(2,:));
	%cost_2=sum(xi1.*xi2*Z.graph.pairwise_str(6,:)  );
	%compute pairwise solution
	v1=[Z.B.pos_inds;Z.B.neg_inds];
	v2=(v1*0)+1;
	y1=Z.B.E([Z.B.pos_inds(:);Z.B.neg_inds(:)],1);
	y2=Z.B.E([Z.B.pos_inds(:);Z.B.neg_inds(:)],2);

	[~,y1]=ismember(y1,Z.B.var_keep);		
	[~,y2]=ismember(y2,Z.B.var_keep);
	v3=gc_sol(y1).*gc_sol(y2);
	gc_pair=sparse(v1,v2,v3,Z.B.N_e,1);
	%compute unary
	v1=Z.B.var_keep;
	v2=(v1*0)+1;
	v3=gc_sol;
	gc_unary=sparse(v1,v2,v3,Z.B.N_d,1);%var_keep(gc_sol);
	%join terms togetehr
	sol=[gc_unary(:);gc_pair];	
	red_cost=inst_cost+(obj'*sol(Z.B.inds_keep));
	cost=inst_cost+(Z.ILP.OBJ'*sol);
	if(0==1)
		[sol_red_1,red_cost_1,flag_1]=intlinprog(obj,Z.ILP.matlab_integ_red,A,B,[],[],lb,ub,Z.ILP.matlab_opt);

		red_cost_1a=inst_cost+red_cost;
        	sol_1a=Z.M.mat_sol_2_full*sol_red_1;
        	cost_1a=inst_cost+(Z.ILP.OBJ'*sol_1a);	
		if(abs(red_cost_1-red_cost)>.0001)
			disp('prob here ')	
			red_cost
			red_cost_1
			save('redProb')
			pause	
		end	
	end

	%cost
	%pause
end
if(inst_cost<0 && sum(sol)<0.5)%check if the solution is empty and have negative instance cost.  THis means thaty ou have an issue
	disp('you have an issue here.  you cant have negative offset and permissible empty solution')
	pause
end
%Z.B.inds_force_on
%disp('done itt')
