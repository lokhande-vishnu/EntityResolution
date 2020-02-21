

%addpath to interface
addpath('../jy_fun');
%




%load path to solver Zeus
addpath('../../../');
num_thresh=20;
%load problem instance

my_times_base=[];
my_times_new=[];
my_bounds_base=[];
my_bounds_new=[];
use_orig=0;
for(i=1:1039)
	rand('twister',0);
	%load(['../examples/multi_person/output/out_thresh_doi_',jy_pad_int(i,4),'_',num2str(num_thresh),'.mat'],['H_comp']);
	load(['../examples/multi_person/output/out_thresh_doi_',jy_pad_int(i,4),'_',num2str(num_thresh),'_use_orig_',num2str(use_orig),'.mat'],['H_comp']);
	my_bounds_base=[my_bounds_base;max(H_comp.H_reg_doi.lb),min(H_comp.H_reg_doi.lp),min(H_comp.H_reg_doi.ub)];
	my_bounds_new=[my_bounds_new;max(H_comp.H_thresh_doi.lb),min(H_comp.H_thresh_doi.lp),min(H_comp.H_thresh_doi.ub)];


	my_times_base=[my_times_base;sum(H_comp.H_reg_doi.timer,1)];
	my_times_new=[my_times_new;sum(H_comp.H_thresh_doi.timer,1)];
	%jy_out_val('sum(H_comp.H_reg_doi.timer)',sum(H_comp.H_reg_doi.timer))
	%jy_out_val('sum(H_comp.H_thresh_doi.timer)',sum(H_comp.H_thresh_doi.timer))
	%jy_out_val('max(H_comp.H_reg_doi.lb),max(H_comp.H_rot_doi.lb)',[max(H_comp.H_reg_doi.lb),max(H_comp.H_thresh_doi.lb)])
	%jy_out_val('max(H_comp.H_reg_doi.lp),max(H_comp.H_rot_doi.lp)',[min(H_comp.H_reg_doi.lp),min(H_comp.H_thresh_doi.lp)])
	%jy_out_val('max(H_comp.H_reg_doi.ub),max(H_comp.H_rot_doi.ub)',[min(H_comp.H_reg_doi.ub),min(H_comp.H_thresh_doi.ub)])
	%save('../output/tst_out.mat','H');
end
([num_thresh,sum(my_times_new(:,2))/sum(my_times_base(:,2))])

[~,my_order]=sort(-my_times_base(:,2));
K=1;
sum(my_times_new(my_order(1:K),2))/sum(my_times_base(my_order(1:K),2))

