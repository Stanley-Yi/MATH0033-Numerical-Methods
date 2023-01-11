%% Setup
close all, clear all, clc
format long, format compact     % Display all digits and condense output in the command window
tol=1e-10;                           % Tolerance for stopping criteria
nmax=50;                           % Maximum number of iterations (in case stopping criterion is not achieved)
fs=18;                              % Setting font size in plots (the default is a bit small)
f=@(x) (x/2)-sin(x)+(pi/6)-(sqrt(3)/2);                 % Define the function f
set(groot,'defaulttextfontsize',fs);
set(groot,'defaultaxesfontsize',fs);
set(groot,'defaultLineLineWidth',2) % Set line width in plots (the default is a bit thin)
set(groot,'defaultContourLineWidth',2) % Set line width in contour plots (the default is a bit thin)
set(0,'DefaultLegendAutoUpdate','off') % Stops Matlab adding extra legend entries automatically
pauses = 0;  % If 1 then pause between plots, if 0 then don't

%% (a) 
% Bisection
x=linspace(-pi,pi,101);               % Define a set of x values for plotting

figure                              % Create a new figure
plot(x,f(x))          % Plot f
grid on
xlabel('x')
ylabel('f(x)')
hold on                             % Hold allows you to overlay things on top of the current plot

[zero,res,niter,itersb]=bisection(f,0,3,tol,nmax)    % Run the bisection method with a=-1, b=1, and tol and nmax as specified above

for i=1:10
    scatter(itersb(i,1),f(itersb(i,1)),'r')   % Plot the first 10 iterates from bisection
    if pauses, pause, end     % pausing after each one if pauses=1
end

%%
% 
%  On the above, I chose interval [0,3] for the intial data.
% 
%% (b) 
% Newton method
df = @(x) (1/2)-cos(x);            % Define the derivative df/dx

x_a=pi;   % Initial guess of \alpha
x_b=-pi/2;   % Initial guess \beta

x=linspace(-pi,pi,101);               % Define a set of x values for plotting
figure
plot(x,f(x),'DisplayName','f')          % Plot f
grid on
xlabel('x')
ylabel('f(x)')
hold on 

[zero_a,res_a,niter_a,itersn_a]=newton(f,df,x_a,tol,nmax)    % Run Newton
[zero_b,res_b,niter_b,itersn_b]=newton(f,df,x_b,tol,nmax)    % Run Newton

scatter(itersn_a(1,1),f(itersn_a(1,1)),100,'pk','DisplayName','x_\alpha')
scatter(itersn_b(i,1),f(itersn_b(i,1)), 100,'DisplayName','x_\beta')
legend('Location','best')

for i=2:5
    scatter(itersn_a(i,1),f(itersn_a(i,1)),100,'pk') 
    if pauses, pause, end     % pausing after each one if pauses=1
end

for i=2:10
    scatter(itersn_b(i,1),f(itersn_b(i,1)), 100,'^r')
    if pauses, pause, end     % pausing after each one if pauses=1
end

df_a = df(zero_a)
df_b = df(zero_b)
%%
%  
%  For $\alpha$, it used 5 iterations to converge, but $\beta$ used 27
%  iterations to converge.
%  The reason is that: According to the Theorem 3.4.1., we can know that if
%  $f$ be twice continuously differentiable, and $f'(x) \neq 0$, then the
%  Newton iteration converges to x quadratically. The df_a is
%  (approximatively) equal to 0, so $\alpha$ converges faster.
%  However, according to the Remark 3.4.2., we can know that if $f'(x) =
%  0$, the convergence of Newton's method is only linear, that is why
%  $\beta$ converges slower.
% 


%% (c) 
% Fixed point iterations
phi=@(x) x-2*((f(x))/(df(x)));               

x=linspace(-pi,pi,101);               % Define a set of x values for plotting
figure
plot(x,f(x),'DisplayName','f')          % Plot f
grid on
xlabel('x')
ylabel('f(x)')
hold on 

[fixp,res,niter,itersfp1]=fixpoint(phi,x_b,tol,nmax) 

for i=1:4
    scatter(itersfp1(i,1),f(itersfp1(i,1)),500,'.r')
    if pauses, pause, end     % pausing after each one if pauses=1
end


%%
%  
%  This modified Newton method used 4 iterations to make $\beta$
%  convergence.
% 
