clc;
clear all;

%load path to LP/ILP cplex solver if I want to use it
%addpath('/data/shaofeiw/ibm/ILOG/CPLEX_Studio128/cplex/matlab/x86-64_linux')
%addpath('/mnt/hdd-4tb/lokhande_files_4T/SOFTWARE/ibm_ilog/cplex/matlab/x86-64_linux')

%load path to Zeus
addpath('../../../jy_fun');
addpath('../../../master');
addpath('../../../rounding');
addpath('../../../pricing');
addpath('../../../body');
addpath('../pre_process')
addpath('../interface_col_gen')
addpath('../../../');%load path to Zeus.m
%name 'input and output files'
%my_file_input='../prob_inst/F.mat';
%my_file_input='../prob_inst/csv_example/F.mat';

my_file_input='./running.mat'
my_file_output= '../../../../music20k/results_2/H_sample_23_running.mat'

load(my_file_input);

[ub,cols_keep]=solve_ilp(T, G);
H.sol=T.X(:,cols_keep);
v1=[];
v2=[];
n_out=size(H.sol,2);
for(i=1:n_out);
    inds_on=find(H.sol(:,i));
    v1a=repmat(inds_on(:),[numel(inds_on),1]);
    v2a=jy_copy_col(inds_on(:),numel(inds_on));
    v1=[v1;v1a];
    v2=[v2;v2a];
end
v3=(v2*0)+1;
H.sol_sf=sparse(v1,v2,v3,G.B.Nd,G.B.Nd);
%computing second alternative format output
H.sol_cluster=H.sol*([1:size(H.sol,2)]');

save(my_file_output,'H');
