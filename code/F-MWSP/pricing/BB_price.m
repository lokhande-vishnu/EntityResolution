function [nX,nTheta,n_aux,n_resid_list,residual,G]=BB_price(m,C,lambda,psi,G,pricing_option)
%
%Purpose
%
%	Apply pricing.  THis is set up to operate under Subset row inequaliteis
%
%
%Input
%
%	m:  sub-problem index
%
%	C: usual
%
%	lambda: usual
%
%	psi:  usual
%
%	G:  usual
%
%	pricing_option:  0 means use heuristic otherwise use exact
%
%
%Output
%
%	nX:  matrix of the new columns
%
%	nTheta:  list of Thete terms assocaited wit hthe new columns
%
%	n_aux:  matrix of the aux info for hte new columns
%
%	n_resid_list:  residuals of hte new columns
%
%	residual:  total residual for hte sub-problem.  THis is the lowest reduced cost UNLESS it is tracking and we generate muliple tracks for the sub-problem that can be consistent.  
%
%	G: updated terms
%

%inital list of columns found
nX=[];
nTheta=[];
n_aux=[];
n_resid_list=[];
residual=0;

heur_list=-inf;%initialize  the heuristic
branch_list=zeros(1,G.B.Nd);%initial branch

init_resid_list=[];%this holds the residual list
best_decode_resid=0;%initialize value of lowest reduced cost solution

step=1;%count steps
while(min(heur_list)<-G.opt.epsilon)%if the heuristic is negative contine 
	%pop this heuristic from the top of the stack
	this_heur=heur_list(1);
	%if this heuritic is worse than the current lowest reduced cost solution then break
	if(this_heur>=best_decode_resid)
		break;
	end

	%load the current branch.  then remove that branch and heuristic from the queue
	this_branch=branch_list(1,:);
	branch_list=branch_list(2:end,:);
	heur_list=heur_list(2:end);
	%grab the column with lowest reduced cost.  This provides a lot of other columns potentially
	[heur,residual,new_x_opt,new_x_list,new_aux_list,new_theta_list,new_resid_list,ind_opt,used_psi,G]=get_column(m,this_branch,C,lambda,psi,G,pricing_option);
	%store the initial list of lower bounds
	if(numel(init_resid_list)==0)
		init_resid_list=new_resid_list;
	end


  	%if the heuristic is an upper bound of the lowest reduced cost column then we dont continue 
	if(heur>best_decode_resid || heur==0)
		continue;
	end
	%get heuristic to be assoiciated with each child branch
	[heur_new]=fix_price(this_branch,heur,C,psi,used_psi);
	%get the  reduced cost of the current produced solution
	[upper_resid]=fix_price(new_x_opt,heur,C,psi,used_psi);
	%get the reduced cost of each produced solution
	new_resid_list=fix_group(new_x_list,new_resid_list,C,psi,used_psi);
	%only retrain indexes with neative reduced cost
	inds_keep=find(new_resid_list<-G.opt.epsilon);
	nX=[nX,new_x_list(:,inds_keep)];
	nTheta=[nTheta;new_theta_list(inds_keep)];
	n_aux=[n_aux,new_aux_list(:,inds_keep)];
	n_resid_list=[n_resid_list;new_resid_list(inds_keep)];
	%grab and hold the optimal solution's value
	if(numel(n_resid_list)>0)
		best_decode_resid=min(n_resid_list);    
	end
	%FIX PRICING DONE
    
	if(   (upper_resid-heur_new)<G.opt.epsilon || numel(psi)<0.5)%if no triples exist then terminate 
		continue;
	else
		save('ddtt')
		disp('I should not be here in hte cell section')
		pause
		%do a split operation
		[new_branch_list,new_heur_list]=do_split(heur_new,new_x_opt,this_branch,C,lambda,psi,G,used_psi);
		%add the new branches that are created
		if(numel(new_heur_list)>0)
            		branch_list=[branch_list;new_branch_list];
			heur_list=[heur_list;new_heur_list];
			[heur_list,new_order]=sort(heur_list);
			heur_list=heur_list(:); 
			branch_list=branch_list(new_order,:);
		end
	end
	%if no heurisitc is in the list then we termninate
	if(numel(heur_list)<0.5)
		break;
	end	
	step=step+1;
end
%residual is the current heuristic
residual=this_heur;

if(G.B.do_sum_list==1)%if we are in the ituation of tracking we add the costs computed
	init_resid_list(init_resid_list<this_heur)=this_heur;%this line sets tracks with cost below truth to the vlaue of hte truth
	residual=sum(init_resid_list.*(double(init_resid_list<=0)));%this selects the tracks with negatrive reduced cost
end

