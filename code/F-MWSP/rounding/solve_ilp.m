function     [ilp_val,cols_keep]=solve_ilp(T,G)

%Purpose
%
%	Produce a rounded solution
%
%Input
%
%	T: Holds current columns
%	G: holds problem structure
%
%Output
%
%	ilp_val:  Objective of integer solution
%	cols_keep:  grabs the columns selected in the ILP solution
%
%Formulate ILP
OBJ=T.Theta;%- is because matlab minimizes not maximizes the value of LPs
A=T.X;
B=ones(size(A,1),1);
%Solve ILP
[gamma,ilp_val]=intlinprog(OBJ,[1:numel(OBJ)],A,B,[],[],0*OBJ,(0*OBJ)+1,G.opt.int_lin_prog_opt );

%grab solution
cols_keep=find(gamma>0.5);

