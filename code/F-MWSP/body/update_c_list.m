function [C,did_add]=update_c_list(G,gamma,T,C)

%Purpose
%
%	determine which subset row inequalities to add.  Add them to C
%
%Input
%
%	G:  usual
%	gamma:     primal solution 
%	T:  usual
%	C: usual
%
%Output
%
%	C: augmented from input
%
%	did_add:  set to one if a new row has been created
%

%Idenfify fractional columns
cols_keep=find(gamma.*(1-gamma)>.001);
cols_keep=cols_keep(:)';
inds_keep=find(sum(T.X(:,cols_keep),2)>0.5);
unique_pairs=[];
%-------------------
%identify unique pairs.  These are all pairs of detections each of which is associated with a single set
for(col=cols_keep)
	%get the active detections
	inds_incl=find(T.X(:,col)==1);
	%get the list of pairs
	tmp1=repmat(inds_incl,[numel(inds_incl),1]);
	tmp2=jy_copy_col(inds_incl,numel(inds_incl));
	%add the detection pairs
	unique_pairs=[unique_pairs;[tmp1,tmp2]];
end
%apply unique term.
unique_pairs=sort(unique_pairs')';
unique_pairs=unique(unique_pairs,'rows');
%--------------------------
%create list of unique triples

%Add all detections to each pair
uni_trip=repmat(unique_pairs,[numel(inds_keep),1]);
tmp=jy_copy_col(inds_keep,size(unique_pairs,1));
uni_trip=[uni_trip,tmp];
uni_trip=sort(uni_trip')';%uni trip
uni_trip=unique(uni_trip,'rows');
%get all removes of trip needing to be removed
remove1=(uni_trip(:,1)==uni_trip(:,2));
remove2=(uni_trip(:,2)==uni_trip(:,3));
remove3=(uni_trip(:,1)==uni_trip(:,3));
uni_trip=uni_trip( ( remove3+ remove1+remove2 ) <0.5,:);
uni_trip=unique(uni_trip,'rows');
%Create matrix to hold triples

v1=jy_copy_col([1:size(uni_trip,1)],3);
v2=reshape(uni_trip',[numel(uni_trip),1]);
v3=(v2*0)+1;
C3_mat_candid=sparse(v1,v2,v3,max(v1),G.B.Nd);
%identify how much each triple is violated.  Retain only violated ones
mat2useA=(C3_mat_candid*T.X>=1.99);
amount_viol=mat2useA*gamma;
inds_retain=find(amount_viol>1.0001);
C3_mat_candid=C3_mat_candid(inds_retain,:);
uni_trip=uni_trip(inds_retain,:);

amount_viol=amount_viol(inds_retain);
%iif there are more triples to add then slots available then add the most violated ones
if(numel(inds_retain)>G.opt.max_trip_add)
    [amount_viol,order]=sort(-amount_viol);
    order=order(1:G.opt.max_trip_add);
    C3_mat_candid=C3_mat_candid(order,:);
    uni_trip=uni_trip(order,:);
end
%add to C the new rows
did_add=0;
if(numel(uni_trip)>0.5)
	did_add=1;

	C.c3_list=[C.c3_list;uni_trip];
	C.c3_mat=[C.c3_mat;C3_mat_candid];

end

