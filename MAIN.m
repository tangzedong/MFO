% This MATLAB R2014b code is for EVOLUTIONARY MULTITASKING across minimization problems. 
% For maximization problems, multiply objective function by -1.

% For suggestions please contact: Abhishek Gupta (Email: abhishekg@ntu.edu.sg or
% agup839@aucklanduni.ac.nz or abhi.nitr2010@gmail.com)

clear
close all
%% Example 1 - (40-D Rastrigin, 30-D Ackley)
% % Rastrigin function definition
n=7;
Tasks(1).dims=n;
M=orth(randn(n,n));
Tasks(1).fnc=@(x)Perm(x,M,zeros(1,1));
Tasks(1).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
Tasks(1).Lb=-n*ones(1,n);
Tasks(1).Ub=n*ones(1,n);
% Ackley function definition
n=7;
Tasks(2).dims=n;
M=orth(randn(n,n));
Tasks(2).fnc=@(x)Ackley(x,M);
Tasks(2).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
Tasks(2).Lb=-50*ones(1,n);
Tasks(2).Ub=50*ones(1,n);

%% Example 2 - (50-D Sphere, 30-D Weierstrass)
% Sphere function definition
% n=30;
% Tasks(1).dims=n;
% M=eye(n,n);
% Tasks(1).fnc=@(x)Rastrigin(x,M);
% Tasks(1).fncLS = @(fnc,x0,option)fminunc(fnc,x0,option);
% Tasks(1).Lb=-50*ones(1,n);
% Tasks(1).Ub=50*ones(1,n);
% % Rastrigin function definition
% n=30;
% Tasks(2).dims=n;
% M=orth(randn(n,n));
% Tasks(2).fnc=@(x)Ackley(x,M);
% Tasks(2).fncLS = @(fnc,x0,option)fminunc(fnc,x0,option);
% Tasks(2).Lb=-50*ones(1,n);
% Tasks(2).Ub=50*ones(1,n);

%% Example 3 - (40-D Rastrigin, 50-D Ackley, 20-D Sphere)
% % Rastrigin function definition
% n=40;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Rastrigin(x,M);
% Tasks(1).Lb=-5*ones(1,n);
% Tasks(1).Ub=5*ones(1,n);
% % Ackley function definition
% n=50;
% Tasks(2).dims=n;
% M=orth(randn(n,n));
% Tasks(2).fnc=@(x)Ackley(x,M);
% Tasks(2).Lb=-32*ones(1,n);
% Tasks(2).Ub=32*ones(1,n);
% % Sphere function definition
% n=20;
% Tasks(3).dims=n;
% M=eye(n,n);
% Tasks(3).fnc=@(x)Sphere(x,M);
% Tasks(3).Lb=-100*ones(1,n);
% Tasks(3).Ub=100*ones(1,n);

%% Example 5 - Mirror Functions - (30-D Rastrigin, 30-D Rastrigin)
% % Rastrigin function definition
% n=30;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Rastrigin(x,M);
% Tasks(1).Lb=-5*ones(1,n);
% Tasks(1).Ub=5*ones(1,n);
% Tasks(2) = Tasks(1);

%% Example 6 - Mirror Functions - (30-D Ackley, 30-D Ackley)
% % Ackley function definition
% n=30;
% Tasks(1).dims=n;
% M=orth(randn(n,n));
% Tasks(1).fnc=@(x)Ackley(x,M);
% Tasks(1).Lb=-32*ones(1,n);
% Tasks(1).Ub=32*ones(1,n);
% Tasks(2)=Tasks(1);

%% Example 7 - Discrete optimization problem
% cvrpdata = CVRPread('A-n65-k9.vrp');
% Tasks(1).dims=cvrpdata.dim - 1;
% Tasks(1).fnc=@(x)CVRP(x,cvrpdata);
% Tasks(1).fncLS=@(fnc, x0, option)CVRPlocalsearch(x0, cvrpdata);
% Tasks(1).Lb=zeros(1,cvrpdata.dim);
% Tasks(1).Ub=ones(1,cvrpdata.dim);
% 
% qapdata = QAPread('CHR22B.txt');
% Tasks(2).dims=qapdata.nnodes;
% Tasks(2).fnc=@(x)QAP(x,qapdata);
% Tasks(2).fncLS=@(fnc, x0, option)QAPlocalsearch(x0, qapdata);
% Tasks(2).Lb=zeros(1,qapdata.nnodes);
% Tasks(2).Ub=ones(1,qapdata.nnodes);

%% Calling the solvers
% For large population sizes, consider using the Parallel Computing Toolbox
% of MATLAB.
% Else, program can be slow.
pop = 100; % population size
gen = 500; % generation count
selection_pressure = 'roulette wheel'; % choose either 'elitist' or 'roulette wheel'
p_il = 1; % probability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp = 0.4; % random mating probability
data_MFEA=MFEA(Tasks,pop,gen,selection_pressure,rmp,p_il);

% task_for_comparison_with_SOO = 2;
% data_SOPSO=SOPSO(Tasks(task_for_comparison_with_SOO),task_for_comparison_with_SOO,pop,gen,selection_pressure,p_il);     
