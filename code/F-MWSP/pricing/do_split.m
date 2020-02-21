function         [new_branch_list,new_heur_list]=do_split(this_heur,this_sol,this_branch,C,lambda,psi,G,used_psi)

%Purpose
%
%	split a branch along a tribple.  This creates many branches.  
%
%
%Input
%
%	this_heur:  lower bound associated with branch all child branches to be created
%	
%	this_sol:  the solution as a single column associated with the heuristic
%
%	this_branch:  description of the branch,  0 for not modeled, 1 for must include,-1 must not include 
%
%	C:  usual
%
%	lambda:  usual
%
%	psi:  usual
%
%	G:  usual
%
%	used_psi:  list of psi terms modeled in the sub-problem
%
%Output
%
%
%	new_branch_list:  holds the new branches
%
%	new_heur_list:  hold the new heuristics
%
%
%

%%%%%%
%Identify the triple associated with the greatest violation
%grab detections in the solution
dets_in_branch=find(this_branch>0.5);
dets_in_sol=find(this_sol>0.5);

psi(used_psi)=0;% sets to zero the psi terms that are modeled in the priceing
q1=(sum(C.C3_mat*(this_branch(:)>0.5),2)>1.5);
q2=(sum(C.C3_mat*(this_sol(:)>0.5),2)>1.5);
[~,split_trip]=max((q2-q1).*psi;% computes the triple that is most violated
%%%%%%

%grab the detections that compose the triple
d1=C.C3_list(split_trip,1);
d2=C.C3_list(split_trip,2);
d3=C.C3_list(split_trip,3);

%Holds the possible ways of splitting 
Y1=[-1,-1, -1,-1,  1, 1,  1,1];
Y2=[-1,-1,  1, 1, -1, -1, 1,1];
Y3=[-1, 1, -1, 1, -1, 1, -1,1];
Y=[Y1',Y2',Y3'];

new_branch_list=[];
new_heur_list=[];
%itterate over possible ways of splitting.  Add them if they are valid
for(i=1:8)
	y=Y(i,:);%  grab way of splitting
	if ( y(1) ~= this_branch(d1)  && this_branch(d1)~=0) %if y(1) is inconsistent with branch Stop ittearion
		continue
	end
	if(  y(2) ~= this_branch(d2)  && this_branch(d2)~=0)%if y(2) is inconstent with branch stop itteraiton
		continue
	end
	if(  y(3) ~= this_branch(d3)  &&  this_branch(d3)~=0)%if y(3) is incocnsitent with branch stop itteration
		continue
	end
	%copy branch from original
	new_branch=this_branch;
	%Set the inclusion/exclusion for the branch
	new_branch(d1)=y(1);
	new_branch(d2)=y(2);
	new_branch(d3)=y(3);
	%add the branch and heuristics to the list
	new_branch_list=[new_branch_list;new_branch];
	new_heur_list=[new_heur_list;this_heur];

end

