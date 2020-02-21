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

my_file_input='../../../../NVCRdataset/results_1/F_sample.mat';
my_file_output='../../../../NVCRdataset/results_1/H_sample.mat';

%load input file to produce F.  
F1=load(my_file_input);%produces F

%set my_opts to empty.  We can alter this if needed
my_opts=[];
my_opts.use_heur=1;%set to 1 if you want to  use heuristic pricing only.  set to 2 if you want heuristic pricing followed by optimal pricing. 0 is only exact pricing

my_opts.num_thresh=0;  %this sets the number of percentiles used in the new DOI formulation.  Set to zero to use the original terms.  
%call pre-process to make problem instance
my_opts.max_trip_add=0;
disp('starting pre-process')
my_opts.do_no_dom_neigh_pricing=1;
rand('twister',0);
%Preprocess start
uni_inds=unique([F1.pairwise(:,[1]);F1.pairwise(:,2);F1.neigh(:)]);
F=F1;
[~,F.pairwise(:,[1,2])]=ismember(F1.pairwise(:,[1,2]),uni_inds);
[~,F.neigh]=ismember(F1.neigh,uni_inds);

%Preprocess

G=pre_process_thresh_low_mem(F,my_opts);

disp('starting Zeus')
G.opt.solve_ilp_prior=0;
G.opt.save_each_itt=0;
%call Zeus to solve optimization
G.opt.max_cols_add=20;
H=Zeus_K_doi_partial(G);
%postprocess

num_nodes=max(double(uni_inds));
sol_new=spalloc(double(num_nodes),size(H.sol,2),sum(H.sol(:)));
sol_new(double(uni_inds),:)=H.sol;
H.sol=sol_new;
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
H.sol_sf=sparse(v1,v2,v3,num_nodes,num_nodes);
%computing second alternative format output
H.sol_cluster=H.sol*([1:size(H.sol,2)]');


%save result
save(my_file_output,'H');
%3249
