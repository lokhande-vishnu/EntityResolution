
function [heur_new,extra_term]=fix_price(x_vec,heur,C,psi,used_psi)

%Purpose
%
%	update the heuristic by including the upper bound component.  
%	This gives the actual reduced cost of the vector
%Input
%
%	x_vec:  column
%	heur:  lower bound on reduced
%	C:  usual
%	psi:  usual
%	use_psi: psi terms already modeled in the pricing problem
%
%Output
%
%
%	heur_new :  Actual Reduced Cost.  heuristic updated with the active psi terms;
%	extra_term= amount added to the heur to get heur_new
%
%

heur_new=heur;%intintalize true reduced cost to heuristic
extra_term=0;
psi(used_psi)=0;%set to zero any psi terms incldued in the heuristic
if(numel(psi)>0)
	extra_term=psi(:)'*(sum(C.c3_mat*(x_vec(:)>0.5),2)>1.5);%grab sum of  active psi terms
	heur_new=heur+extra_term;  %add in the extra term to get the  actual reduced cost
end
