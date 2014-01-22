function [f,g] = myfun(x)
%  x : initial point
%  f : function value
%  g : gradient
%
% example,
%  [f,g] = myfun([1 1 1 1 1 1]);
%
f = 3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 +x(2)+ x(3)^2 + x(4)^2 -x(4)+ x(5)^4+6*x(6)^2+3*x(6);
g = zeros(length(x),1);
if nargout > 1
   g(1) = 6*x(1)+2*x(2);
   g(2) = 2*x(1)+2*x(2)+1;
   g(3) = 2*x(3);
   g(4) = 2*x(4)-1;
   g(5) = 4*x(5)^3;
   g(6) = 12*x(6)+3;
end
