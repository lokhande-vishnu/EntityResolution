
function T2=call_reprice(G,candid_col);

%PURPOSE
%
%	Call the reprice operation on T2
%
%INPUT
%
%
%	G:  usual	
%
%	candid_col:  binary vector of size N_d
%
%OUTPUT
%
%
%	T2:  same form as T
%	
%
%write command
out_portion='T2=';
in_portion='(G,candid_col);';
fun_portion=G.B.reprice_command;

my_command=[out_portion,fun_portion,in_portion];


%call evaluation code
eval(my_command);

%keep only those with negative cost
if(numel(T2.Theta)>0)
	new_inds=find(T2.Theta<0);
	
	T2.X=T2.X(:,new_inds);
	T2.aux_info=T2.aux_info(:,new_inds);
	T2.Theta=T2.Theta(:,new_inds);
	
end

