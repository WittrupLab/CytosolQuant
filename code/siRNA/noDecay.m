function f = noDecay(x, theta)
%%
%   function with expontial decay
%
%   f(x) =  I(x>\theta_1-\theta_3)I(x < \theta_1)  +
%           I(x > \theta_1) * (  \theta_2)
%   x     - (n x 1) data positions
%   theta - (6 x 1) parameters
%

%%
x = x(:);
theta(1) = exp(theta(1))-2;
theta(2)  =exp(theta(2));
theta(3) = exp(theta(3))/(1+exp(theta(3)));
f = zeros(length(x),1);
x_d = x - theta(1);
index = x_d>0;
f(index) =  theta(2);
index = (x > theta(1) - theta(3) ).* (x < theta(1));
if(sum(index)>0)
   f(index==1) =   f(index==1) +  theta(2) * (x(index==1) - theta(1) + theta(3));
end
