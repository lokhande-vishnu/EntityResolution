

%addpath to interface
addpath('../jy_fun');
%




%load path to solver Zeus
addpath('../../../');
num_rot=10;
%load problem instance

my_times_base=[];
my_times_new=[];
my_bounds_base=[];
my_bounds_new=[];
for(i=1:3844)
	rand('twister',0);
	load(['../examples/multi_person/output/out_rot_doi_',jy_pad_int(i,4),'_',num2str(num_rot),'.mat'],['H_comp']);

	my_bounds_base=[my_bounds_base;max(H_comp.H_reg_doi.lb),min(H_comp.H_reg_doi.lp),min(H_comp.H_reg_doi.ub)];
	my_bounds_new=[my_bounds_new;max(H_comp.H_rot_doi.lb),min(H_comp.H_rot_doi.lp),min(H_comp.H_rot_doi.ub)];


	my_times_base=[my_times_base;sum(H_comp.H_reg_doi.timer,1)];
	my_times_new=[my_times_new;sum(H_comp.H_rot_doi.timer,1)];
	jy_out_val('sum(H_comp.H_reg_doi.timer)',sum(H_comp.H_reg_doi.timer))
	jy_out_val('sum(H_comp.H_rot_doi.timer)',sum(H_comp.H_rot_doi.timer))
	jy_out_val('max(H_comp.H_reg_doi.lb),max(H_comp.H_rot_doi.lb)',[max(H_comp.H_reg_doi.lb),max(H_comp.H_rot_doi.lb)])
	jy_out_val('max(H_comp.H_reg_doi.lp),max(H_comp.H_rot_doi.lp)',[min(H_comp.H_reg_doi.lp),min(H_comp.H_rot_doi.lp)])
	jy_out_val('max(H_comp.H_reg_doi.ub),max(H_comp.H_rot_doi.ub)',[min(H_comp.H_reg_doi.ub),min(H_comp.H_rot_doi.ub)])
	%save('../output/tst_out.mat','H');
end
