

function [max_val,max_trip,DOI]=get_DOI_col(G,C,this_col,DOI)

%Purpose
%
%	Execute the command to get DOI for a given column
%
%Input
%
%	G:  usual
%	C:  usual
%	this_col:  one column of size G.B.Nd,1
%	this_aux_info:  one column's aux info
%
%Output
%
%	max_val:  upper bound on dual values associated with a given col
%		size Nd,1
%	max_trip:
%		size N_trip,1
%

%combine the terms into command
out_portion='[max_val,max_trip,DOI]=';
in_portion='(G,C,this_col,DOI);';
fun_portion=G.B.DOI_command;
my_command=[out_portion,fun_portion,in_portion];


%execute command
eval(my_command);



