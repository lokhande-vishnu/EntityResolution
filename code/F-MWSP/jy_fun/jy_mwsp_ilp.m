function [inds_integ,timer_ILP]=jy_mwsp_ilp(Z)

%Purpose
%
%	return indexes that should be integral.  DOne using ILP solver
%
%Input
%
%	Z
%
%
%Output
%
%
%	inds_integ= output integs
%
%


timer_ILP=0;
M_worry=[];
inds_integ=[];%initalize inds_integ to empty
if( Z.B.N_pos >0)

	v1=repmat([1:Z.B.N_pos]',[2,1]);%jy_copy_col([1:size(E_worry,1)]',2);
	v2=Z.B.E(Z.B.pos_inds,:);
	v2=v2(:);
	v3=(v2*0)+1;
	M_worry=sparse(v1,v2,v3,max(v1),Z.B.N_d);
end
%add in the hard cliques
M_worry=[M_worry;Z.B.hard_cliques];
if(Z.B.N_pos>0 || numel(Z.B.hard_cliques)>0)
	M_worry(:,Z.B.inds_force_on)=0;

	A=M_worry(:,Z.B.var_keep);
	
	A=A(1.5<  sum(A,2)  ,: );
	B=(A(:,1)*0)+1;
	C=-ones(numel(Z.B.var_keep),1);
	lb=0*C;
	ub=1+lb;
	[~,tmp]=ismember(Z.B.inds_force_on,Z.B.var_keep);
	ub(tmp)=0;
	my_opt=optimoptions('intlinprog','Display','off');
	timer_ILP=tic();
	[sol,val]=intlinprog(C,[1:numel(C)],A,B,[],[],lb,ub,my_opt);
	timer_ILP=toc(timer_ILP);
	sol(tmp)=1;
	%sol(Z.B.inds_force_on)=1;
	inds_integ=find(sol<0.5);
	inds_integ=Z.B.var_keep(inds_integ);
end
