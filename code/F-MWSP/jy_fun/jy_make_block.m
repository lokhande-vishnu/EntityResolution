function S1=jy_make_block(N)
S=[0;1];
for(n=1:N-1)
    x0=zeros(size(S,1),1);
    x1=ones(size(S,1),1);
    S=[S,x0;S,x1];
end

if(N<1)
    S=[];
end


