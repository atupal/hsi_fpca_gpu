%  ----------------------  test without weights  ----------------------

n = 21;

tvec = linspace(0,1,n)';

xvec = sin(2*pi*tvec);

sigerr = 0.1;

yvec = xvec + randn(n,1)*sigerr;

argvals    = tvec;
y          = xvec;
wtvec      = [];
covariates = [];
fdnames    = [];

nbasis = n+2;
basisobj = create_bspline_basis([0,1],nbasis);
fdParobj = basisobj;

lambda = 1e-4;

fdParobj = fdPar(basisobj, 2, lambda);

% tic;
% [fdobj, df, gcv, beta, SSE] = ...
%                  smooth_basis(argvals, y, fdParobj);
% toc;
tic;
[fdobj, df, gcv, beta, SSE] = ...
                 smooth_basis_old(argvals, y, fdParobj);
toc;

%  test with covariate

zvec = zeros(n,1);
zvec(11) = 1;
beta = 0.5;

yvec = yvec + zvec*beta;
y = yvec;

lambda = 1e-4;
fdParobj = fdPar(basisobj, 2, lambda);

[fdobj, df, gcv, beta, SSE] = ...
                 smooth_basis(argvals, y, fdParobj);
[fdobj, df, gcv, beta, SSE] = ...
                 smooth_basis(argvals, y, fdParobj, ...
                              'covariates', zvec);
[fdobj, df, gcv, beta, SSE] = ...
                 smooth_basis(argvals, y, fdParobj, ...
                              'covariates', zvec, 'method', 'qr');

beta

figure(1)
plotfit_fd(y, argvals, fdobj)

%  test with multiple records

xmat = [sin(2*pi*tvec),cos(2*pi*tvec)];

sigerr = 0.1;

ymat = xmat + randn(n,2)*sigerr;

beta = [-1,1];
ymat(:,1) = ymat(:,1) + zvec*beta(1);
ymat(:,2) = ymat(:,2) + zvec*beta(2);

y = ymat;

lambda = 1e-4;
fdParobj = fdPar(basisobj, 2, lambda);

[fdobj, df, gcv, beta, SSE] = ...
                 smooth_basis(argvals, y, fdParobj, ...
                                 'covariates', zvec);

beta

plotfit_fd(y, argvals, fdobj)

%  ----------------------  test with weight vector  ----------------------

N = 50;  %  number of curves
n = 21;  %  number of sampling points

%  range of argument values and argument vector

T = 20;
argvals = linspace(0,T,n)';

%  set up variances for residuals

%  set up a random residual variance distribution

% sigvar = 0.4;
% resvar = exp(randn(n,1).*sigvar);

%  set up a smooth residual variance distribution

resvar = exp(sin(2*pi*argvals/T));

%  the weight vector is the reciprocal of the variance vector

wtvec0 = 1./resvar;
wtvec  = wtvec0;

figure(1)
subplot(2,1,1)
plot(argvals, resvar, 'o', [0,T], [1,1], 'r:')
subplot(2,1,2)
plot(argvals, wtvec0, 'o', [0,T], [1,1], 'r:')

%  a saturated order 4 spline basis for the data

nbasis = n + 2;
basisobj = create_bspline_basis([0,T], nbasis);

%  random coefficients for the smooth functions

sigcoef = 1;
coefmat = randn(nbasis,N)*sigcoef;

xfd = fd(coefmat, basisobj);

xmat = eval_fd(argvals, xfd);

%  define error matrix and noisy data

sigerr = 0.4;
errmat = mvnrnd(zeros(n,1), diag(resvar), N)'.*sigerr;

res0 = errmat;

ymat = xmat + errmat;
y = ymat;

figure(1)
subplot(2,1,1)
plot(argvals, errmat)
subplot(2,1,2)
plot(argvals, ymat)

%  loop through log10 smoothing parameter values, 
%  compuing GCV, degrees of freedom, and root mean squared error for each
%  value

lnlam    = -3:0.25:0;
gcvsave  = zeros(length(lnlam),1);
dfsave   = gcvsave;
RMSEsave = gcvsave;
for i=1:length(lnlam)
  fdPari = fdPar(basisobj, 2, 10^lnlam(i));
  %  smoothing without weighting
%   [yfdi, dfi, gcvi] = smooth_basis(argvals, ymat, fdPari);
  %  smoothing with weighting
  [yfdi, dfi, gcvi] = smooth_basis(argvals, ymat, fdPari, 'weight', wtvec0);
  gcvsave(i)  = sum(gcvi);
  dfsave(i)   = dfi;
  resi        = eval_fd(argvals, yfdi) - xmat;
  RMSEsave(i) = sqrt(mean(resi(:).^2));
end

%  plot the results

figure(2)
subplot(3,1,1)
phdl = plot(lnlam, gcvsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} GCV(\lambda)')
subplot(3,1,2)
phdl = plot(lnlam, RMSEsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} RMSE(\lambda)')
subplot(3,1,3)
phdl = plot(lnlam, dfsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}(\lambda)')
ylabel('\fontsize{13} df(\lambda)')

%  smooth for a specific log10 smoothing parameter

lnlam = -0.5;
fdParobj = fdPar(basisobj, 2, 10^lnlam);
[yfd1, df1, gcv1] = ...
    smooth_basis(argvals, ymat, fdParobj);
[yfd1, df1, gcv1] = ...
    smooth_basis(argvals, ymat, fdParobj, ...
                 'weight', wtvec0, 'method', 'qr');
% [yfd1, df1, gcv1] = ...
%     smooth_basis(argvals, ymat, fdParobj, wtvec0);

figure(3)
plotfit_fd(ymat, argvals, yfd1)

res1 = ymat - eval_fd(argvals, yfd1);

wtvec1 = mean(res1.^2,2);

figure(4)
plot(sigerr.*resvar, wtvec1, 'o')

wtfdParobj = fdPar(basisobj, 2, 1);
wtfd1 = smooth_basis(argvals, wtvec1, wtfdParobj);
wtvec1 = eval_fd(argvals, wtfd1);

figure(4)
plot(sigerr.*resvar, wtvec1, 'o')

RMSE1 = sqrt(mean(res1(:).^2));

disp([RMSE1, df1, sum(gcv1)])

lnlam    = -3:0.25:0;
gcvsave  = zeros(length(lnlam),1);
dfsave   = gcvsave;
RMSEsave = gcvsave;
for i=1:length(lnlam)
  fdPari = fdPar(basisobj, 2, 10^lnlam(i));
  [yfdi, dfi, gcvi] = ...
      smooth_basis(argvals, ymat, fdPari, 'weight', wtvec1);
  gcvsave(i)  = sum(gcvi);
  dfsave(i)   = dfi;
  resi        = eval_fd(argvals, yfdi) - xmat;
  RMSEsave(i) = sqrt(mean(resi(:).^2));
end

%  plot the results

figure(3)
subplot(3,1,1)
phdl = plot(lnlam, gcvsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} GCV(\lambda)')
subplot(3,1,2)
phdl = plot(lnlam, RMSEsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} RMSE(\lambda)')
subplot(3,1,3)
phdl = plot(lnlam, dfsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}(\lambda)')
ylabel('\fontsize{13} df(\lambda)')

lnlam = -2;
fdParobj = fdPar(basisobj, 2, 10^lnlam);
[yfd2, df2, gcv2] = ...
    smooth_basis(argvals, ymat, fdParobj, 'weight', wtvec1);

figure(4)
plotfit_fd(ymat, argvals, yfd2)

res2 = eval_fd(argvals, yfd2) - xmat;

RMSE2 = sqrt(mean(res2(:).^2));

disp([RMSE2, df2, sum(gcv2)])

%  ----------------------  test with weight matrix  ----------------------

N = 50;  %  number of curves
n = 21;  %  number of sampling points

%  range of argument values and argument vector

T = 20;
argvals = linspace(0,T,n)';

%  set up covariance matrix for residuals

addpath('../fdVariance')

alpha = 1;
beta  = 0.02;
nu    = 4;

sigma = maternSig(beta,nu,n,T,alpha);

figure(1)
surf(sigma)

wtmat = inv(sigma);
wtmat0 = wtmat./max(max(wtmat));

figure(2)
surf(wtmat0)

%  convert weight vector to weight matrix

wtmat = diag(wtvec0);

nbasis = n + 2;
basisobj = create_bspline_basis([0,T], nbasis);

sigcoef = 1;
coefmat = randn(nbasis,N)*sigcoef;

xfd = fd(coefmat, basisobj);

xmat = eval_fd(argvals, xfd);

sigerr = 0.4;
errmat = mvnrnd(zeros(n,1), sigma, N)'.*sigerr;

ymat = xmat + errmat;
y = ymat;

figure(2)
subplot(2,1,1)
plot(argvals, errmat)
subplot(2,1,2)
plot(argvals, ymat)

lnlam    = -3:0.25:0;
gcvsave  = zeros(length(lnlam),1);
dfsave   = gcvsave;
RMSEsave = gcvsave;
for i=1:length(lnlam)
  fdPari = fdPar(basisobj, 2, 10^lnlam(i));
%   [yfdi, dfi, gcvi] = smooth_basis(argvals, ymat, fdPari);
  [yfdi, dfi, gcvi] = smooth_basis(argvals, ymat, fdPari, 'weight', wtmat);
  gcvsave(i)  = sum(gcvi);
  dfsave(i)   = dfi;
  resi        = eval_fd(argvals, yfdi) - xmat;
  RMSEsave(i) = sqrt(mean(resi(:).^2));
end

%  plot the results

figure(3)
subplot(3,1,1)
phdl = plot(lnlam, gcvsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} GCV(\lambda)')
subplot(3,1,2)
phdl = plot(lnlam, RMSEsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} RMSE(\lambda)')
subplot(3,1,3)
phdl = plot(lnlam, dfsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}(\lambda)')
ylabel('\fontsize{13} df(\lambda)')

lnlam = -1.25;
fdParobj = fdPar(basisobj, 2, 10^lnlam);
[yfd1, df1, gcv1] = smooth_basis(argvals, ymat, fdParobj);

res1 = eval_fd(argvals, yfd1) - xmat;

RMSE1 = sqrt(mean(res1(:).^2));

disp([RMSE1, df1, sum(gcv1)])

lnlam    = -3:0.25:0;
gcvsave  = zeros(length(lnlam),1);
dfsave   = gcvsave;
RMSEsave = gcvsave;
for i=1:length(lnlam)
  fdPari = fdPar(basisobj, 2, 10^lnlam(i));
  [yfdi, dfi, gcvi] = ...
      smooth_basis(argvals, ymat, fdPari, 'weight', wtmat0);
  gcvsave(i)  = sum(gcvi);
  dfsave(i)   = dfi;
  resi        = eval_fd(argvals, yfdi) - xmat;
  RMSEsave(i) = sqrt(mean(resi(:).^2));
end

%  plot the results

figure(4)
subplot(3,1,1)
phdl = plot(lnlam, gcvsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} GCV(\lambda)')
subplot(3,1,2)
phdl = plot(lnlam, RMSEsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('')
ylabel('\fontsize{13} RMSE(\lambda)')
subplot(3,1,3)
phdl = plot(lnlam, dfsave, '-o');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{13} log_{10}(\lambda)')
ylabel('\fontsize{13} df(\lambda)')

lnlam = -12;
fdParobj = fdPar(basisobj, 2, 10^lnlam);
[yfd2, df2, gcv2] = smooth_basis(argvals, ymat, fdParobj, 'weight', wtmat);

res2 = eval_fd(argvals, yfd2) - xmat;

RMSE2 = sqrt(mean(res2(:).^2));

disp([RMSE2, df2, sum(gcv2)])

