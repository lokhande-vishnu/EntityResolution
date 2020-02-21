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

my_file_input='../../../../music20k/results_2/F_sample.mat';
my_file_output='../../../../music20k/results_2/H_sample_22.mat';

%load input file to produce F.  
F=load(my_file_input);%produces F

%set my_opts to empty.  We can alter this if needed
my_opts=[];
my_opts.use_heur=1;%set to 1 if you want to  use heuristic pricing only.  set to 2 if you want heuristic pricing followed by optimal pricing. 0 is only exact pricing
                  
my_opts.num_thresh=3;  %this sets the number of percentiles used in the new DOI formulation.  Set to zero to use the original terms.  


%call pre-process to make problem instance
my_opts.max_trip_add=0;
disp('starting pre-process')
%my_opts.max_sz_node=5; % uncomment for affiliations dataset
%my_opts.do_no_dom_neigh_pricing=1; % uncomment for affiliations dataset
rand('twister',0);
G=pre_process_thresh_low_mem(F,my_opts);

disp('starting Zeus')
G.opt.solve_ilp_prior=0;
G.opt.save_each_itt=1;
G.save_running_name='running.mat'
%call Zeus to solve optimization
%G.opt.max_cols_add=20; % uncomment for affiliations dataset
H=Zeus_K_doi_partial(G);
%save result
save(my_file_output,'H');

