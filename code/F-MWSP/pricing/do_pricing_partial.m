function     [G,T,residual,inds_new,max_time_pricing,current_index]=do_pricing_partial(G,C,T,lambda,psi,pricing_option,current_index)%Purpose
%
%	Iterate over pricing problems and solve each.  	
%
%Input
%
%
%	G: usual
%
%	C: usual
%
%	T: usual
%
%	psi: usual
%
%	pricing_option:  0 means we use heuristic and otherwise use exact solver
%
%	current_index
%Output
%	
%	G:  This may be modfied to change the solver info
%	
%	T:  usual with added cols
%
%	residual: sum of all residuals in all pricing problems.  Produces a valid lower boudn when added to LP
%
%	inds_new:  indexes of new cols
%
%	max_time_pricing:  maximum time spent on a pricing problem
%
%	current_sub-problem being explored
%
%Important Variables
%
%	new_T:  holds all the columns
%		
%	resid_list:  holds the residuals 
%

%initialize the new columns
new_T=[];
new_T.Theta=[];
new_T.X=[];
new_T.aux_info=[];

%initilalize residual list
residual=0;
resid_list=[];

%
max_time_pricing=0;
num_add=0;
%my_order=randperm(G.B.num_pricing);
count=1;
my_order=[current_index:G.B.num_pricing,1:(current_index-1)];
for(i=my_order)

	current_index=i+1;
	m=G.B.sub_prob_order(i);
%for(m=1:G.B.num_pricing)
	%solve the sub-problelem
	%jy_out_val('size(resid_list,1)',size(resid_list,1))
	this_time_pricing=tic();
    	[nX,nTheta,n_aux,resid_list_part,residual_this,G]=BB_price(m,C,lambda,psi,G,pricing_option);
	max_time_pricing=max([toc(this_time_pricing),max_time_pricing]);
	
	%if residual is negative then add the columns
	if(residual_this<-G.opt.epsilon)
        	new_T.X=[new_T.X,nX];
       		new_T.Theta=[new_T.Theta;nTheta];
        	new_T.aux_info=[new_T.aux_info,n_aux];
        	resid_list=[resid_list;resid_list_part];% add residual info
      		if(size(resid_list,1)>G.opt.max_cols_add)% &&  sum(resid_list)<-100)
			disp('breaking due to adding enough cols')	
			break;
		end	      
    	end 
	%add the residual to lower bound 
	residual=residual+residual_this;
	count=count+1;
end
if(current_index>G.B.num_pricing)
	current_index=1;
end
jy_out_val('count sub prob check',count)
%get theunique cols
[~,rows_keep]=unique(new_T.X','rows');
new_T.Theta=new_T.Theta(rows_keep);
new_T.X=new_T.X(:,rows_keep);
new_T.aux_info=new_T.aux_info(:,rows_keep);
resid_list=resid_list(rows_keep);

%if there are a too many cols generated then only keep the most violated.  
if(numel(resid_list)>G.opt.max_cols_add)
	[~,my_order]=sort(resid_list);
	rows_keep=my_order(1:G.opt.max_cols_add);
	new_T.Theta=new_T.Theta(rows_keep);
	new_T.X=new_T.X(:,rows_keep);
	new_T.aux_info=new_T.aux_info(:,rows_keep);
	resid_list=resid_list(rows_keep);
end
%identify the new indexes 
inds_new=numel(T.Theta)+(1:numel(new_T.Theta));
%Add to T the new_T
T.Theta=[T.Theta;new_T.Theta];
T.X=[T.X,new_T.X];
T.aux_info=[T.aux_info,new_T.aux_info];
