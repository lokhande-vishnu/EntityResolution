clear
%load path to LP/ILP cplex solver if I want to use it
addpath('/mnt/hdd-4tb/lokhande_files_4T/SOFTWARE/ibm_ilog/cplex/matlab/x86-64_linux')
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
%my_file_input='/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/result_11/F_sample_r1.mat';
%my_file_input='/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/result_13/F_sample.mat';
%my_file_input='/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/music20k_sota/result_1/F_sample.mat';
%my_file_input='/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/result_13/F_sample.mat';%'/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/music20k_sota/F_sample.mat';
%my_file_input='../prob_inst/F_too_long_1.mat';%csv_example_drjulian/Full.mat';
%my_file_input='../prob_inst/csv_example_drjulian/affiliations_sample/result_1/F_sample.mat';
%my_file_input='/home/julian/Zeus/examples/entity_rez_vishnu/prob_inst/csv_example_drjulian/result_14/F_sample.mat';


my_file_input='../../../../music20k/results_2/F_sample.mat';
my_file_output='../../../../music20k/results_2/H_sample_21.mat';


% my_file_output='../output/ex_output.mat';
%load input file to produce F.  
F=load(my_file_input);%produces F

%set my_opts to empty.  We can alter this if needed
my_opts=[];
my_opts.use_heur=0;%set to 1 if you want to  use heuristic pricing only.  set to 2 if you want heuristic pricing followed by optimal pricing. 0 is only exact pricing

my_opts.num_thresh=0;  %this sets the number of percentiles used in the new DOI formulation.  Set to zero to use the original terms.  
%call pre-process to make problem instance
my_opts.max_trip_add=0;
disp('starting pre-process')
my_opts.max_sz_node=5;
G=pre_process_thresh_bound(F,my_opts);

disp('starting Zeus')
%call Zeus to solve optimization
H=Zeus_K_doi(G);
%save result
save(my_file_output,'H');
%3249
