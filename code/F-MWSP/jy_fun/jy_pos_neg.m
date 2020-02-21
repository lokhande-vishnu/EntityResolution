function [y1,y2]=jy_pos_neg(x)

y1=x.*(double(x>0));
y2=x.*(double(x<0));