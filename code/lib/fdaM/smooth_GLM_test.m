%  tests for function smooth_basis_LS

%  ----------------  normal link with no covariates  --------------------

n       = 101;
argvals = linspace(0,2*pi,n)';
y0      = sin(2*argvals);
y02     = cos(2*argvals);
sig     = 0.2;
y       = y0 + sig.*randn(n,1);
y       = [y, y02 + sig.*randn(n,1)];

basisobj = create_bspline_basis([0,2*pi],n+2);

basismat = eval_basis(argvals, basisobj);

Lfdobj = int2Lfd(2);
penmat = eval_penalty(basisobj,Lfdobj);

lambda  = 1e-1;
lamRmat = lambda.*penmat;

[coef,Deviance] = glm_fda(basismat, y, 'normal', lamRmat);

fdobj = fd(coef,basisobj);

plotfit_fd(y, argvals, fdobj)

fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev] = ...
            smooth_GLM(argvals, y, fdParobj, 'family', 'normal');
                  
plotfit_fd(y, argvals, fdobj)

df

gcv

SSE

dev

beta

%  ----------------  normal link with covariates  --------------------

n       = 101;
argvals = linspace(0,2*pi,n)';
y0      = sin(2*argvals);
sig     = 0.2;
y       = y0 + sig.*randn(n,1);
sigcov  = 0.1;
covariates = randn(n,1);
beta    = 1;
y       = y + covariates*beta;

basisobj = create_bspline_basis([0,2*pi],11);

basismat = eval_basis(argvals, basisobj);

Lfdobj = int2Lfd(2);
penmat = eval_penalty(basisobj,Lfdobj);

lambda  = 1;
lamRmat = lambda.*penmat;

fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev] = ...
            smooth_GLM(argvals, y, fdParobj, ...
                       'family', 'normal', 'covariates', covariates);
                  
plotfit_fd(y, argvals, fdobj)

beta

df

gcv

SSE

dev

%  ----------------  binomial link with no covariate  --------------------

n       = 501;
argvals = linspace(0,1,n)';
y0      = sin(4*pi*argvals);
sig     = 0.5;
y = zeros(n,1);
y(y0 + sig.*randn(n,1) >= 0.0) = 1;

basisobj = create_bspline_basis([0,1],13);

basismat = eval_basis(argvals, basisobj);

Lfdobj = int2Lfd(2);
penmat = eval_penalty(basisobj,Lfdobj);

lambda  = 1e-4;
lamRmat = lambda.*penmat;

[coef,Deviance] = glm_fda(basismat, y, 'binomial', lamRmat);

fdobj = fd(coef,basisobj);

plotfit_fd(y, argvals, fdobj)

fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev] = ...
                      smooth_GLM(argvals, y, fdParobj, ...
                                       'family', 'binomial');
                  
plot(fdobj)

beta

df

gcv

SSE

dev

stats{:}

%  ----------------  poisson link with no covariates  --------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);
y = exp(y);

basisobj = create_bspline_basis([0,2*pi],53);

lambda = 1e-1;
Lfdobj = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
            smooth_basis_GLM(argvals, y, fdParobj, ...
                             'family', 'poisson');

plotfit_fd(log(y), argvals, fdobj)

beta

df

gcv

SSE

dev

stats{:}


