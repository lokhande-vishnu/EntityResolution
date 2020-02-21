function my_num=jy_pad_int(x,size_min)

%purpose
%
%	take an input integer and covert it to a string with zeros in front so that there is a minimum size of int
%
%input
%
%	x:  positive integer
%	size_min:  minimum size of string
%
%output
%
%	my_num:  string augmented with zeros with original number in front
%
my_num=num2str(x);
n_zeros_add=max([0,size_min-numel(my_num)]);
for(j=1:n_zeros_add)
	my_num=['0',my_num];
end

