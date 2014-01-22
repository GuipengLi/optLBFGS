% Here is an simple example of how to use optLBFGS to solve optimization(min) problem. 
% we need to know the function value and gradient , as myfun return
x0 = ones(6,1);
%
maxIter = 100;
% memsize <= length(x0)
% % memsize : 5~20 will be fine
memsize = 3; 
tic
[X,F,k] =optLBFGS(@myfun,x0,maxIter,memsize)
toc

% another example ,the 2D Rosenbrock "banana" function
% http://www.di.ens.fr/~mschmidt/Software/minFunc.html
D = 20;
x1 = zeros(D,1);
[X1,F1,k1] =optLBFGS(@rosenbrock,x1,200,5);
[X2,F2]  = minFunc(@rosenbrock,x1);
