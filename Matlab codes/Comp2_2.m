%% Computational Problem 2

%% Setup
close all, clear all, clc
format long, format compact
fs = 16;
set(0,'defaulttextfontsize',fs);
set(0,'defaultaxesfontsize',fs);


%% EXERCISE 2


%% (a)
% 

list_N = [5, 10, 20, 40, 80];
rhoBJ = zeros(size(length(list_N)));
rhoBGS = zeros(size(length(list_N)));

for i = 1:length(list_N)
    h = 1 / list_N(i);
    A = (2/h^2)*diag(ones(list_N(i)-1,1)) - (1/h^2)*diag(ones(list_N(i)-2,1),1) - (1/h^2)*diag(ones(list_N(i)-2,1),-1);
    b = transpose(sin(pi*h*(1:list_N(i)-1)));

    D = diag(diag(A));
    L = tril(A) - D;
    U = triu(A) - D;

    B_J = -(D^(-1)) * (L + U);
    B_GS = -(D + L)^(-1) * U;

    rhoBJ(i) = max(abs(eig(B_J)));  % spectral radius is the maximum modulus of eigenvalues
    rhoBGS(i) = max(abs(eig(B_GS)));
end

rhoBJ, rhoBGS


%%
% The two iterative methods will convergent for these values of N, because
% each of spectral radius less than 1. And the Gauss-Seidel method converge
% faster, because every value in rhoBGS less than rhoBJ.
%
% As size N increases, both $\rho\left(B_{\mathrm{J}}\right)$ and
% $\rho\left(B_{\mathrm{GS}}\right)$ are grow up. Thus, I expect that the
% performance of the Jacobi and Gauss-Seidel methods will decline with the
% size N continues to increase.
%

figure

Nvec = [5, 10, 20, 40, 80];
loglog(Nvec, 1 - rhoBJ);
hold;
loglog(Nvec, 1 - rhoBGS);

legend({'rhoBJ', 'rhoBGS'}, 'Location', 'Best');


% get relationship between size N and spectral radius

% for $\rho B_J$
syms a b;
f1 = a*list_N(1) + b == 1 - rhoBJ(1);
f2 = a*list_N(5) + b == 1 - rhoBJ(5);
ab_J = solve(f1, f2, a, b);
vpa(ab_J.a), vpa(ab_J.b)

% for $\rho B_{GS}$
syms a b;
f1 = a*list_N(1) + b == 1 - rhoBGS(1);
f2 = a*list_N(5) + b == 1 - rhoBGS(5);
ab_GS = solve(f1, f2, a, b);
vpa(ab_GS.a), vpa(ab_GS.b)

f_J = @(x) 0.0025 * x + 0.796;
f_GS = @(x) 0.0046 * x + 0.632;

% relationship between size N and spectral radius for both methods

x=linspace(5, 100, 100);               % Define a set of x values for plotting
figure                              % Create a new figure
plot(x, f_J(x))          % Plot f
hold
plot(x, f_GS(x))          % Plot f

legend({'rhoBJ', 'rhoBGS'}, 'Location', 'Best');
grid on
xlabel('x')
ylabel('f(x)')

%%
%
% According to the graph above, we can observe that the function between
% size N and (1 - spectral radius) is a straight line. So, we can deduce
% that their relationship is linear, and the function is: (1 - sr) = aN + b,
% where sr is spectral radius, a and b are parameters.
%
% For $\rho\left(B_{\mathrm{J}}\right)$, we can calculate that a = -0.0025,
% b = 0.204, so the final relationship is sr = 0.0025N - 0.204 + 1 =>
% sr = 0.0025N + 0.796.
%
% For $\rho\left(B_{\mathrm{GS}}\right)$, we can calculate that a =
% -0.0046, b = 0.368, so the final relationship is sr = 0.0046N - 0.368 + 1
% => sr = 0.0046N + 0.632.
%
%
% According to the functions above, as size N increases, the spectral
% radius also increases, so that the number of iterations to achieve a
% fixed solution accuracy will increases as well. Additionally, the
% spectral radius may greater than 1 when N approaches infinity, so the two
% methods will not converge, and the iterations number will also approaches
% infinity.
% 
% The cost of a direct method is generally O(n^3), where n is the size of
% matrix; and for iterative method, the cost of each iteration in general
% is O(n^2). So, the Jacobi and Gauss-Seidel methods cheaper than direct
% method if the number of iterations to be much less than N. However, when
% N is large, iteration number will larger than N if we use Jacobi and
% Gauss-Seidel methods to solve this problem. Therefore, in terms of
% practicality, using the Jacobi and Gauss-Seidel methods for this problem
% is inferior to using direct method.
% 



%% (b)
figure
grid on
hold on

tol = 1e-10;
nmax = 10^5;

list_N = [5, 10, 20, 40, 80];

% Plot resulting solutions of Gauss-Seidel method

for i = 1:length(list_N)
    h = 1 / list_N(i);
    A = (2/h^2)*diag(ones(list_N(i)-1,1)) - (1/h^2)*diag(ones(list_N(i)-2,1),1) - (1/h^2)*diag(ones(list_N(i)-2,1),-1);
    b = transpose(sin(pi*h*(1:list_N(i)-1)));
    x0 = transpose(zeros(1, list_N(i)-1));

    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');

    result = zeros(1, list_N(i) + 1);
    result(2:list_N(i)) = x;
    range = linspace(0, 1, list_N(i) + 1);

    plot(range, result,'LineWidth', 2);
end

% Plot exact solution y(x) on the finest mesh (N = 80)
N = 80;
h = 1 / N;
y = @(x) pi^(-2) * sin(pi * x);
range = linspace(0, 1, N + 1);
result = zeros(1, N + 1);

for i = 2:N
    result(i) = y(i * h);
end

plot(range, result,'LineWidth', 2);

legend({'N=5', 'N=10', 'N=20', 'N=40', 'N=80', 'y(x)'}, 'Location', 'Best');

%%
% According to the graph, as N increases, the solution appears closer to
% the line of exact solution, and the error also becomes smaller. Thus, the
% solutions appear to converge as N increases.
% 

list_N = [5, 10, 20, 40, 80];
error_vect = zeros(1, length(list_N));

y = @(x) pi^(-2) * sin(pi * x);

for i = 1:length(list_N)
    h = 1 / list_N(i);
    A = (2/h^2)*diag(ones(list_N(i)-1,1)) - (1/h^2)*diag(ones(list_N(i)-2,1),1) - (1/h^2)*diag(ones(list_N(i)-2,1),-1);
    b = transpose(sin(pi*h*(1:list_N(i)-1)));
    x0 = transpose(zeros(1, list_N(i)-1));

    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');
    
    % calculate y(x_n)
    y_n = transpose(zeros(1, list_N(i)-1)); 
    for j = 1:list_N(i)-1
        y_n(j) = y(j * h);
    end
    
    error_vect(i) = max(abs(x - y_n));

end

error_vect

figure
hvect = 1./list_N;
loglog(hvect, error_vect)  % e(N)
hold


syms c
cond1 = c * hvect(1)^2 >= error_vect(1);
cond2 = c * hvect(5)^2 >= error_vect(5);
conds = [cond1 cond2];

C = solve(conds, c);
vpa(C)

loglog(hvect, C*hvect.^2);  % E(h)

legend({'e(N)', 'E(h)'}, 'Location', 'Best');


%%
% 
% According to the assumption of question, $e(N) \leq Ch^2$ as h = 1/N $\rightarrow$ 0,
% so we can know that error e(N) is proportional to the step size h = 1/N
% from the graph.
% To verify the assumption, we set p=2, and find C by the graph. Thus,
% C=1.0833396945217224072166573023424. Finally, we can see that the line of
% e(N) lower than the line of E(h), which means $e(N) \leq Ch^2$.
% Therefore, the numerical results support this theoretical estimate.
% 


%% (c)
figure
grid on
hold on

tol = 1e-10;
nmax = 10^5;

list_N = [5, 10, 20, 40, 80];

% Plot resulting solutions of Gauss-Seidel method

for i = 1:length(list_N)
    h = 1 / list_N(i);
    A = (2/h^2)*diag(ones(list_N(i)-1,1)) - (1/h^2)*diag(ones(list_N(i)-2,1),1) - (1/h^2)*diag(ones(list_N(i)-2,1),-1);
    b = transpose(ones(1, list_N(i)-1));  % modify
    x0 = transpose(zeros(1, list_N(i)-1));

    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');

    result = zeros(1, list_N(i) + 1);
    result(2:list_N(i)) = x;
    range = linspace(0, 1, list_N(i) + 1);

    plot(range, result,'LineWidth', 2);
end

% Plot exact solution y(x) on the finest mesh (N = 80)
N = 80;
h = 1 / N;
y = @(x) (1 / 2) * x * (1 - x);
range = linspace(0, 1, N + 1);
result = zeros(1, N + 1);

for i = 2:N
    result(i) = y(i * h);
end

plot(range, result,'LineWidth', 2);

legend({'N=5', 'N=10', 'N=20', 'N=40', 'N=80', 'y(x)'}, 'Location', 'Best');


% Corresponding errors

Nvec = [5, 10, 20, 40, 80];
error_vect = zeros(1, length(Nvec));

y = @(x) (1 / 2) * x * (1 - x);

for i = 1:length(Nvec)
    h = 1 / Nvec(i);
    A = (2/h^2)*diag(ones(Nvec(i)-1,1)) - (1/h^2)*diag(ones(Nvec(i)-2,1),1) - (1/h^2)*diag(ones(Nvec(i)-2,1),-1);
    b = transpose(ones(1, Nvec(i)-1));  % modify
    x0 = transpose(zeros(1, Nvec(i)-1));

    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');
    
    % calculate y(x_n)
    y_n = transpose(zeros(1, Nvec(i)-1)); 
    for j = 1:Nvec(i)-1
        y_n(j) = y(j * h);
    end
    
    error_vect(i) = max(abs(x - y_n));

end

figure
hvect = 1./Nvec;
loglog(hvect, error_vect)


%%
% For this case, there is no longer a relationship of $e(N) \leq Ch^2$ as h
% = 1/N $\rightarrow$ 0. Since we replaced $sin(\pi x)$ with 1, there is no x exists in
% $d^2/dx^2$, so constant C is not proportional to $\max _{x \in[0,1]}\left|y^{\prime \prime \prime \prime}(x)\right|$.
% 