
function [rows_use,rows_use_trip]=select_rows_use(G,T,C,DOI)
%Purpose
%
%	produce rows
%
%output
%
%
%
%	rows_use
%	
%		each row has d,threshold,gain

rows_use=[];
rows_use_trip=[];
n_thresh_default=1+G.DOI_opt.num_percentiles;
for(d=1:G.B.Nd)
	
	%min_val
	inds_keep=find(DOI.max_val_mat(d,:)>abs(G.opt.epsilon));	
	list_vals=DOI.max_val_mat(d,inds_keep);
	%list_vals
	%jy_out_val('min and max vals',[min(list_vals),max(list_vals)])
	%pause
	list_vals=unique(floor(list_vals*10000)/10000);
	if(numel(list_vals)>0.5)
		
		inds_cut=[];
		if(numel(list_vals)<=G.DOI_opt.num_percentiles)
			inds_cut=list_vals;
		else
			if(0.5<G.DOI_opt.num_percentiles)
			
				n_vals=numel(list_vals);
				list_vals=sort(list_vals);
				inc=n_vals/G.DOI_opt.num_percentiles;
				my_keep=1:inc:numel(list_vals);
				my_keep=ceil(my_keep);
				my_keep=[1,my_keep,numel(list_vals)];
				my_keep=unique(my_keep);
				inds_cut=list_vals(my_keep);
				%inds_cut
				%pause
			else
				inds_cut=max(list_vals);
			end
		end
		c1=[];
		c2=[];
		c3=[];
		if(0.5<G.DOI_opt.num_percentiles)
			n_cut=numel(inds_cut);
			c1=d*ones(n_cut,1);
			%c2=inds_cut*.99999;
			c3=.001+(1.00001*(inds_cut-[0,inds_cut(1:end-1)]));
			c2=[inds_cut(1)/2,(.00001+[inds_cut(1:end-1)])];
		else
			c1=d;
			c2=min(list_vals)*.9999;
			c3=(max(list_vals)*1.0001)+.0001;
		end
		rows_use=[rows_use;c1,c2(:),c3(:)];	
		%if(d==76)
		%	inds_cut	
		%	jy_out_val('[c1,c2(:),c3(:)]',[c1,c2(:),c3(:)])
		%end
		%rows_use
		%pause
	end
end

num_trip=size(C.c3_list,1);
if(num_trip>0.5)
	disp('not setup yet')
	pause
	for(c=1:num_trip)
	
	end
end
