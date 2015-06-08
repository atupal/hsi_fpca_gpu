%  test the use of dropind

%  set up a basis

rangeval = [0,1];
nbasis   = 5;
norder   = 4;
dropind  = 1;
breaks   = linspace(0,1,3);
basisobj0 = create_bspline_basis(rangeval, nbasis, norder, breaks);
basisobj  = create_bspline_basis(rangeval, nbasis, norder, breaks, dropind);
% basisobj = create_fourier_basis(rangeval, nbasis, 1, dropind);
% print the basis
display(basisobj);
% plot the basis (this tests a number of other functions, too)
plot(basisobj);

inprod_basis(basisobj0, basisobj0)
inprod_basis(basisobj, basisobj)

%  evaluate the penalty matrix
full(eval_penalty(basisobj,0))
%  set up a fd object
fdobj = fd(randn(4,1),basisobj);
plot(fdobj)

x = linspace(0,1,21)';
y = sin(2*pi*x) + 0.2.*randn(21,1);

ysmth = smooth_basis(x, y, basisobj);

fdParobj = fdPar(basisobj, 2, 1e-4);

ysmth = smooth_basis(x, y, fdParobj);

plot(ysmth)

x = linspace(0,1,21)';
y(:,1,1) = sin(2*pi*x) + 0.5.*randn(21,1);
y(:,2,1) = cos(2*pi*x) + 0.5.*randn(21,1);
y(:,1,2) = sin(2*pi*x) + 0.5.*randn(21,1);
y(:,2,2) = cos(2*pi*x) + 0.5.*randn(21,1);

fdobjx = smooth_basis(x, y(:,:,1), fdParobj);
fdobjy = smooth_basis(x, y(:,:,2), fdParobj);
ncan = 2;
ccafdParx = fdParobj;
ccafdPary = fdParobj;

ccaStruct = cca_fd(fdobjx, fdobjy, ncan, ccafdParx, ccafdPary, 0);
ccaStruct = cca_fd(fdobjx, fdobjy, ncan);

ccaStruct
ccaStruct.corrs
plot(ccaStruct.wtfdx)
plot(ccaStruct.wtfdy)

pcastr = pca_fd(fdobjx, 2, fdParobj, 0);

pcastr = pca_fd(fdobjx, 2);

xfd = fd(randn(4,21),basisobj);

betafd = fd(randn(5,1),basisobj0);

efd = fd(randn(4,21).*0.1,basisobj);

%  this gives an order 12 spline with 13 basis functions for yfd

yfd = betafd.*xfd + efd;

%  this gives an order 4 spline with 5 basis functions for yfd

tfine = linspace(0,1,201)';
xmat = eval_fd(tfine, xfd);
bmat = eval_fd(tfine, betafd);
emat = eval_fd(tfine, efd);
ymat = (bmat*ones(1,21)).*xmat + emat;
yfd = smooth_basis(tfine, ymat, basisobj0);

yfdPar = fdPar(yfd);
size(getcoef(getfd(yfdPar)))

xfdcell{1} = xfd;
betacell{1} = betafd;

fRegressStruct = fRegress(yfd, xfdcell, betacell);

%  plot estimated beta

betahat   = fRegressStruct.betahat;
plotbeta(betahat)
betahatfd = getfd(betahat{1});
bhatmat   = eval_fd(tfine, betahatfd);
plot(tfine, bhatmat, '-', tfine, bmat, 'b--')

yhatfd = fRegressStruct.yhat;
yhatmat = eval_fd(tfine, yhatfd);
for i=1:21
    plot(tfine, yhatmat(:,i), '-', tfine, ymat(:,i), 'b--')
    pause
end


