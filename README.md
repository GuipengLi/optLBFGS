optLBFGS
========

Matlab code for the Limited-memory BFGS (Broyden–Fletcher–Goldfarb–Shanno) algorithm.

1. Limited-memory BFGS (L-BFGS) is an optimization algorithm in the family of quasi-Newton methods that approximates the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm using a limited amount of computer memory. It is a popular algorithm for parameter estimation in machine learning.
http://en.wikipedia.org/wiki/Limited-memory_BFGS

2. I use line search algorithm satisfying strong Wolfe conditions.
More details can be found from Algorithms 3.2 on page 59 in Numerical Optimization, by Nocedal and Wright
http://sentientdesigns.net/math/mathbooks/Number%20theory/Numerical%20Optimization%20-%20J.%20Nocedal,%20S.%20Wright.pdf

3. optLBFGS has similar performance with minFunc( with limited-memory BFGS - default)
