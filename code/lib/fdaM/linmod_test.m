%  linmod test

addpath('../..')

sbasis  = create_fourier_basis([0,1],3);
tbasis0 = create_fourier_basis([0,1],5);
tbasis1 = create_fourier_basis([0,1],7);
ybasis  = create_fourier_basis([0,1],9);
xbasis  = create_fourier_basis([0,1],11);

coef0  = randn(5,1);
stcoef = randn(3,7);

N = 20;
sig = 1;
xcoef = randn(11,N);
ecoef = randn( 9,N).*sig;
xfd   = fd(xcoef,xbasis);
afd   = fd(coef0, tbasis0);
stmat = inprod(sbasis, xfd);
betaxintfd = fd(stcoef'*stmat,tbasis1);

plot(betaxintfd)

yfd = betaxintfd;
yfd = afd + betaxintfd;
efd = fd(ecoef,ybasis);
yfd = afd + betaxintfd + efd;

betacell{1} = fdPar(afd);
betafd = bifd(stcoef, sbasis, tbasis1);
betacell{2} = bifdPar(betafd);

linmodstr = linmod(xfd, yfd, betacell);

alphahatfd = linmodstr.alpha;
betahatfd  = linmodstr.beta;
yhatfd     = linmodstr.yhat;

figure(1)
plot(yfd)
figure(2)
plot(yfd - yhatfd)
 
figure(1)
plot(afd)
figure(2)
plot(alphahatfd)

sfine = linspace(0,1,51)';
tfine = sfine;
betamat    = eval_bifd(sfine, tfine, betafd);
betahatmat = eval_bifd(sfine, tfine, betahatfd);

figure(1)
surf(sfine, tfine, betamat');
figure(2)
surf(sfine, tfine, betahatmat');
 