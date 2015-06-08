
%  test function glm_fda.m

distr = 'normal';
distr = 'binomial';
distr = 'poisson';
distr = 'gamma';
distr = 'inverse gaussian';

nbasis = 13;
basisobj = create_bspline_basis([0,1],nbasis);
Rmat = eval_penalty(basisobj,2);
lambda = 1e-4;
lamRmat = lambda.*Rmat;

N       = 101;
Tvec    = linspace(0,1,N)';
Wtvec   = ones(N,1);
addterm = [];
Xmat = eval_basis(Tvec,basisobj);

Bvec0 = [];

switch distr
case 'normal'
    Ymat = zeros(N,2);
    eta1  = sin(2*pi*Tvec);
    Ymat(:,1) = eta1 + randn(N,1)*0.2;
    eta2  = cos(2*pi*Tvec);
    Ymat(:,2) = eta2 + randn(N,1)*0.2;
case 'binomial'
    M = 5;
    eta1 = 3*sin(2*pi*Tvec)*ones(1,M);
    eta2 = 3*cos(2*pi*Tvec)*ones(1,M);
    Umat1 = rand(N,M);
    Umat2 = rand(N,M);
    Ymattmp = zeros(N,2);
    for i=1:N
        Uveci = Umat1(i,:);
        Ymati = 1./(1+exp(-eta1(i,:)));
        Ycnti = zeros(1,M);
        Ycnti(Ymati >= Uveci) = 1;
        Ymattmp(i,1) = sum(Ycnti);
        Uveci = Umat2(i,:);
        Ymati = 1./(1+exp(-eta2(i,:)));
        Ycnti = zeros(1,M);
        Ycnti(Ymati >= Uveci) = 1;
        Ymattmp(i,2) = sum(Ycnti);
    end
    Ymat = cell(1,2);
    Ymat{1} = Ymattmp;
    Ymat{2} = M*ones(N,2);
case 'poisson'
    eta   = sin(2*pi*Tvec);
    Ymat = exp(eta+ randn(N,1)*0.2);
case 'gamma'
    eta   = sin(2*pi*Tvec).^2 + 1;
    Ymat = 1./(eta + randn(N,1)*0.2);
case 'inverse gaussian'
    eta   = sin(2*pi*Tvec).^2 + 1;
    Ymat = 1./sqrt(eta + randn(N,1)*0.2);
otherwise
    error('Distribution name is invalid.');
end

[Bvec,Deviance] = ...
              glm_fda(Xmat, Ymat, distr, lamRmat, Wtvec, Bvec0, addterm);

disp(['Deviance = ',num2str(Deviance)])
Yhat = Xmat*Bvec;

figure(1)
subplot(1,1,1)
switch distr
    case 'normal'
        plot(Tvec, Ymat, 'o', Tvec, Yhat, '-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'binomial'
        Phat = 1./(1+exp(-Yhat));
        plot(Tvec, Ymat{1}./Ymat{2}, 'o', Tvec, Phat, 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
   case 'poisson'
        plot(Tvec, Ymat, 'o', Tvec, exp(Yhat), 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'gamma'
        plot(Tvec, Ymat, 'o', Tvec, 1./Yhat, 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'inverse gaussian'
        plot(Tvec, Ymat, 'o', Tvec, 1./sqrt(Yhat), 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
end

%       'dfe'       degrees of freedom for error
%       's'         theoretical or estimated dispersion parameter
%       'sfit'      estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals
%       'residp'    Pearson residuals
%       'residd'    deviance residuals
%       'resida'    Anscombe residuals

disp(['Deviance = ',num2str(Deviance)])

disp(['degrees of freedom for error    = ',num2str(stats.dfe)])
disp(['dispersion parameter            = ',num2str(stats.s)])

