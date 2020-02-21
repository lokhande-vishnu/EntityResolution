

function T2=reprice_obj(G,x);

%purpose
%
%	reprice column
%
%input
%
%
%	G:  usual
%		
%	x:  a column
%
%output
%
%	T2: same as T format
%

%greedily remove elements that make the objective worse until nothing makes the objective worse or only one element reamains
while(sum(x)>1.5)
	
	tmp=x.*(G.Q.unary+(sum(G.Q.aff*x,2)));
	[val,ind]=max(tmp);
	%remove element if needed
	if(val>0.00001)
		x(ind)=0;
	else
		%nothing to remove so break 
		break
	end
	
end
%calclulate cost
this_theta=(0.5*sum(sum(G.Q.aff(x>0.5,x>0.5))))+sum(G.Q.unary(x>0.5))+G.Q.inst_cost;

T2=[];
if(this_theta<0)%add new elements if there is the cost is negative.  
	T2.X=x;
	T2.aux_info=[0;0];
	T2.Theta=this_theta;
else			%add nothing
	T2.Theta=[];
	T2.X=[];
	T.aux_info=[];
end
