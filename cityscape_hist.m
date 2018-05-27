function [new_n,new_x]=cityscape_hist(n,x)

new_x=nan(1,length(x)*2-1);
j=1;
x_step=mode(diff(x));
new_n=nan(1,length(x)*2-1);
for i=1:length(x)-1
    new_x(j)=x(i);
    new_n(j)=n(i);
    j=j+1;
    new_x(j)=x(i)+x_step;
    new_n(j)=n(i);
    j=j+1;
end