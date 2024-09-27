function [ x ] = invlogit(x)

x = 1./(1+exp(-x));

end
