function output=jy_copy_col(input_vec,num_copies)

input_vec=reshape(input_vec,[numel(input_vec),1]);
input_vec=repmat(input_vec,[1,num_copies]);
input_vec=input_vec';
output=reshape(input_vec,[numel(input_vec),1]);