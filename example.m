% Here is an simple example of how to use optLBFGS to solve optimization(min) problem. 
% we need to know the function value and gradient , as myfun return
x0 = ones(6,1);
%
maxIter = 100;
% memsize <= length(x0)
% % memsize : 5~20 will be fine
memsize = 3; 
tic
[X,F,k] =optLBFGS(@myfun,x0,maxIter,memsize);
toc
