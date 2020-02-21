

function LP=form_master_LP_DOI_K(G,T,C,DOI,rows_use)

if(C.c3_mat>0.5)
	disp('not set up for this yet')
	pause
end

LP.obj=[T.Theta(:);rows_use(:,3)];
LP.lb=0*LP.obj;
LP.ub=LP.lb+inf;

num_con=size(rows_use,1);
LP.B=ones(num_con,1);

v1=[];
v2=[];
v3=[];

for(r=1:num_con)

	d=rows_use(r,1);
	v2a=find(DOI.max_val_mat(d,:)>rows_use(r,2));
	v1a=r+(0*v2a);
	v3a=1+(0*v2a);
	v1=[v1;v1a(:)];
	v2=[v2;v2a(:)];
	v3=[v3;v3a(:)];
end

A1=sparse(v1,v2,v3);
A2=-speye(num_con);
LP.A=[A1,A2];
%disp('svaingin in inner dont do')
%save('badCHeck')
%disp('paused')
%pause
if(1<0.5)
	min_red_terms=T.Theta+A1'*rows_use(:,3);
	disp('dododo')
	if(min(min_red_terms)<0)
		disp('bad HERE')
		save('badRed')
		pause
	end
	
end
