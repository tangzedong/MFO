testfun(1).dims = n;
testfun(1).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(1).fnc=@(x)Sumdiff(x,M);
testfun(1).Lb=-50*ones(1,n);
testfun(1).Ub=50*ones(1,n);

testfun(2).dims = n;
testfun(2).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(2).fnc=@(x)Schwefel2(x,M);
testfun(2).Lb=-50*ones(1,n);
testfun(2).Ub=50*ones(1,n);

testfun(3).dims = n;
testfun(3).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(3).fnc=@(x)Rosenbrock(x,M);
testfun(3).Lb=-50*ones(1,n);
testfun(3).Ub=50*ones(1,n);

testfun(4).dims = n;
testfun(4).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(4).fnc=@(x)Perm(x,M);
testfun(4).Lb=-50*ones(1,n);
testfun(4).Ub=50*ones(1,n);

testfun(5).dims = n;
testfun(5).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(5).fnc=@(x)Griewangk(x,M);
testfun(5).Lb=-50*ones(1,n);
testfun(5).Ub=50*ones(1,n);

testfun(6).dims = n;
testfun(6).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(6).fnc=@(x)Axis(x,M);
testfun(6).Lb=-50*ones(1,n);
testfun(6).Ub=50*ones(1,n);

testfun(7).dims = n;
testfun(7).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(7).fnc=@(x)Alpine(x,M);
testfun(7).Lb=-50*ones(1,n);
testfun(7).Ub=50*ones(1,n);

testfun(8).dims = n;
testfun(8).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(8).fnc=@(x)Rastrigin(x,M);
testfun(8).Lb=-50*ones(1,n);
testfun(8).Ub=50*ones(1,n);

testfun(9).dims = n;
testfun(9).fncLS=@(fnc, x0, option)fminunc(fnc, x0, option);
testfun(9).fnc=@(x)Ackley(x,M);
testfun(9).Lb=-50*ones(1,n);
testfun(9).Ub=50*ones(1,n);