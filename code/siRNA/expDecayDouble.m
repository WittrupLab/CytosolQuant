function f = expDecayDouble(x, theta)
%%
%   function with expontial decay with extra jump
%
%   f(x) = I(x>\theta_1-\theta_6)I(x < \theta_1)  +
%          I(x > \theta_1) * (  exp( (\theta_2 -  x^{\theta_5})/\theta_3) + \theta_4)
%          + I(x> \theta_1 + \theta_6) *  exp( (\theta_7 -  x^{\theta_9})/\theta_8)
%   x     - (n x 1) data positions
%   theta - (10 x 1) 
%                   (1) start position
% D: 2019-02-14
%%
x = x(:);
theta(1) = exp(theta(1))-2;
theta(4) = exp(theta(4));
theta(5) = 1+3.*exp(theta(5))/(1+exp(theta(5)));

theta(3) = 4/(-log(0.5))^(1/theta(5)) + exp(theta(3));
f = zeros(length(x),1);
x_d = x - theta(1);
index = x_d>0;
f(index) = exp( theta(2)  - (x_d(index)/theta(3)).^theta(5)) + theta(4);

theta(6) = exp(theta(6))/(1+exp(theta(6)));
index = (x > theta(1) - theta(6) ).* (x < theta(1));
if(sum(index)>0)
   f(index==1) =   f(index==1) + (exp( theta(2) ) + theta(4)) * (x(index==1) - theta(1) + theta(6));
end


theta(7) = exp(theta(7)) + theta(1);
theta(9) = exp(theta(9));
theta(10) = 1+3.*exp(theta(10))/(1+exp(theta(10)));
theta(9)  = 4/(-log(0.5))^(1/theta(10)) + exp(theta(9));
x_d = x - theta(7);
index = x_d>0;
if(sum(index)>0)
    f(index) = f(index) + exp( theta(8)  - (x_d(index)/theta(9)).^theta(10)); 
end

