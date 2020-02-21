

function  my_mat=construct_matrix_rotation_person(G,DOI,r)


Ne=size(G.Q.E,1);
v1=[G.Q.E(:,1);G.Q.E(:,2)];
v2=[1:Ne,1:Ne]';
v3=-repmat(G.Q.E(:,3),[2,1]);
%v3b=G.Q.E(:,3);


do_keep=v3*0;
do_keep(v3>0)=1;%all attractive indexes are kept

left_is_neck=ismember([G.Q.E(:,1);G.Q.E(:,2)],G.Q.F.B.dets_of_fixed);
right_is_neck=ismember([G.Q.E(:,2);G.Q.E(:,1)],G.Q.F.B.dets_of_fixed);

do_keep(left_is_neck)=1;
do_keep(right_is_neck)=1;

inv_list=DOI.inv_rotation_list(r,:);

j1=[G.Q.E(:,1);G.Q.E(:,2)];
j2=[G.Q.E(:,2);G.Q.E(:,1)];
do_keep(inv_list(j1)<inv_list(j2))=1;

v3=v3.*do_keep;


my_mat=sparse(v1,v2,v3,G.B.Nd,Ne);

