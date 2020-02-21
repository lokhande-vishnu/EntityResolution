

function [max_val,max_trip,DOI]=get_DOI_col_rot(G,C,this_col,DOI,r)

%Purpose
%
%	Execute the command to get DOI for a given column
%
%Input
%
%	G:  usual
%	C:  usual
%	this_col:  one column of size G.B.Nd,1
%	DOI:  
%	r:  ordering being manipulated
%
%%

%Output
%
%	max_val:  upper bound on dual values associated with a given col
%		size Nd,1
%	max_trip:
%		size N_trip,1
%	DOI.  updated as needed
%

%combine the terms into command
out_portion='[max_val,max_trip,DOI]=';
in_portion='(G,C,this_col,DOI,r);';
fun_portion=G.B.DOI_command;
my_command=[out_portion,fun_portion,in_portion];


%execute command
eval(my_command);



