function inds_integ=jy_get_integ_greedy(Z)

%purpose
%
%	use a greedy proceadure to construct integral indexes
%
%Input Z
%
%	Z:  standard ILP form
%
%Output
%
%	inds_integ:  indexes that should be integral
%

%first grab any indexes taht are in the hard cliques

%M_worry will hold all of  the relavant constraints

%grab negative indexes
M_worry=[];
inds_integ=[];%initalize inds_integ to empty
if( Z.B.N_pos >0)
	
	v1=repmat([1:Z.B.N_pos]',[2,1]);%jy_copy_col([1:size(E_worry,1)]',2);
	v2=Z.B.E(Z.B.pos_inds,:);
	v2=v2(:);
	v3=(v2*0)+1;
	M_worry=sparse(v1,v2,v3,max(v1),Z.B.N_d);
	if(max(M_worry(:))>1.1)
		disp('error hree')
		pause
	end
end
%jy_out_val('[max(M_worry)]',max(M_worry(:)))
%add in the hard cliques
M_worry=[M_worry;Z.B.hard_cliques];

row_under=[];
if(numel(M_worry)>0.5)%remove the rows containing only only a forced on entry
	row_under=(M_worry(:,1)*0)+1;

	if(numel(Z.B.inds_force_on)>0.5)
		%rows_keep=find(0.5>sum(  M_worry(:,Z.B.inds_force_on) ,2));
		rows_remove=find(     0.5<sum(  M_worry(:,Z.B.inds_force_on) ,2));
		M_worry(rows_remove,:)=0;
		%jy_out_val('rows_remove_1',rows_remove)
		rows_under(rows_remove)=0;
	end
	
end
if(numel(M_worry)>0.5)%remove rows taht contain only one entry
	
	rows_remove=find(1.5 >sum(M_worry,2)) ;
	M_worry(rows_remove,:)=0;
	row_under(rows_remove)=0;
	%jy_out_val('rows_remove2', rows_remove)
end
while(numel(M_worry)>0.5)%itterate until no problematic indexes exist
	
	%[max_val,m]=max(sum(M_worry(row_under>0,:),1));%grab modei
	[max_val,m]=max(sum(M_worry,1));
	if(max_val<0.5)
		break
	end
	inds_integ=[inds_integ,m];%add m to integral list
	M_worry(:,m)=0;
	%rows_poss=find(rows_under>0.5);
	%tmp=find(1.5 >sum(M_worry(rows_poss,:),2));
	%rows_remove=rows_poss(tmp);
	%row_under(rows_remove)=0;
	rows_remove=find(1.5 >sum(M_worry,2)) ;
	M_worry(rows_remove,:)=0;
	%jy_out_val('[m,max_val]',[m,max_val])
end
%jy_out_val('numel(Z.B.var_keep),numel(inds_integ)', [numel(Z.B.var_keep),numel(inds_integ)])
if(numel(inds_integ)==numel(Z.B.var_keep))
	disp('ducks')
	save('duclk')
	pause;
end
%save('ddd')
%jy_out_val('(Z.B.inds_force_on',Z.B.inds_force_on)
%disp('paused')
%pause
