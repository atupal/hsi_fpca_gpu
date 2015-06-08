%  test FEM functions

%  ----  Example 1:  Boundary a right triangle, no interior points  -------

p = [0 0 1; 1 0 0]';

e = [1 2 3; 2 3 1; 0 0 0; 1 1 1; 1 2 3; 1 1 1; 0 0 0]';

t = [1; 2; 3]';

%  create the basis object

basisobj = create_FEM_basis(p, e, t);

plot(basisobj)

params = getbasispar(basisobj);

evalarg.pts = p';
evalarg.z   = [0 1 0]';

u = FEM(evalarg, params);

pdeplot(p, e, t, 'xydata', u, 'colormap', 'hot')

penmat = full(FEMpen(basisobj, int2Lfd(0)))

penmat = full(FEMpen(basisobj, int2Lfd(1)))

penmat = full(eval_penalty(basisobj, int2Lfd(1)))

fdobj = fd(u, basisobj);


%  Example 2:  A right triangle subdivided into four right triangles

p = [0   0  1.0  0    0.5  0.5; ...
     1.0 0  0    0.5  0    0.5];

e = [1.0    2.0    3.0    4.0    5.0    6.0; ...
     4.0    5.0    6.0    2.0    3.0    1.0; ...
     0      0      0      0.5    0.5    0.5; ...
     0.5    0.5    0.5    1.0    1.0    1.0; ...
     1.0    2.0    3.0    1.0    2.0    3.0; ...
     1.0    1.0    1.0    1.0    1.0    1.0; ...
     0      0      0      0      0      0];

t = [3     1     2     4; ...
     6     4     5     5; ...
     5     6     4     6; ...
     1     1     1     1];
 
basisobj = create_FEM_basis(p, e, t);

plot(basisobj)

params = getbasispar(basisobj);

evalarg.pts = p';
evalarg.z   = [0 0 0 1 1 1]';

u0 = FEM(evalarg, params);

pdeplot(p, e, t, 'xydata', u0, 'colormap', 'hot')

penmat = full(FEMpen(basisobj, int2Lfd(0)))

penmat = full(FEMpen(basisobj, int2Lfd(1)))

fdobj = fd(u, basisobj);

plot(fdobj)

fdParobj = fdPar(basisobj, 1, 1e-4);

argvals = p';
y = u0;

fdobj = smooth_FEM_basis(argvals, y, fdParobj);

plot(fdobj)

%  Note:  supplying smooth_FEM_basis with noisy data to be interpolated
%  on to points is risky ... in this example the delaunay triangulation
%  wasn't the same as t, and the results were wrong.
%  Perhaps it would be better to keep this process of data interpolation
%  outside of the smoothing function so that it can be checked.

%  In this version noisy data is supplied as a noisy coefficient vector u0
%  The equations are set up by assema.

lambda = 1e-2;  % a heavy level of smoothing
lambda = 1e-5;  % a light level of smoothing

tctrvec = pdeintrp(p, t, u0);

[K, M, F] = assema(p, t, lambda, 1.0, tctrvec);

u = (K + M)\F;

fdobj = fd(u,basisobj);

%  call smooth_FEM_fd

fdobj = fd(u0, basisobj);
fdsmthobj = smooth_FEM_fd(fdobj, lambda);
plot(fdsmthobj)
xlabel('\fontsize{13} X')
ylabel('\fontsize{13} Y')
title(['\fontsize{13} lambda = ',num2str(lambda)])

%  ------------------  single element  ---------------------

p = [0,0; 2,0; 1,sqrt(3)];
t = [1, 2, 3];
e = [1 2 3; 2 3 1; 0 0 0; 1 1 1; 1 2 3; 1 1 1; 0 0 0]';

[nodes, nodemesh] = makenodes(p,t);

basisobj = create_FEM_basis(p, e, t);

coef = ones(6,1);
fdobj = fd(coef,basisobj);

ngrid = 5;
ngrid = 21;
X = linspace(0,2,ngrid)';
Y = linspace(0,sqrt(3),ngrid)';

plot_FEM(fdobj,X,Y)

u = coef;
tn = X;
al2 = Y;

[uxy,tn,al2,al3] = tri2grid(p,t,u,tn,al2,al3);

evalmat = eval_FEM_fd(X,Y,fdobj);

plot3(X, Y, evalmat, 'o')

data = zeros(6,2);
data(:,1) = (1:6)';
data(:,2) = randn(6,1);

lambda = 1e-4;

newfdobj = smooth_FEM_fd(data,fdobj,lambda);

Xmat = X*ones(1,ngrid);
Ymat = ones(ngrid,1)*Y';
Xvec = Xmat(:);
Yvec = Ymat(:);

evalmat = eval_FEM_fd(Xvec, Yvec, newfdobj);

evalmat = reshape(evalmat,ngrid,ngrid);

surf(evalmat)

plot_FEM(newfdobj, X, Y)

plot(newfdobj)

%  ------------------  single element with additional points  ------------

p = nodes;

p = [p; [1, sqrt(3)/2]];

t = delaunay(p(:,1), p(:,2));

triplot(t,p(:,1),p(:,2))

pdemesh(p',e',t')

nt = size(t,1);
np = size(p,1);

[nodes, nodemesh] = makenodes(p,t);
nNodes = size(nodes,1);

basisobj = create_FEM_basis(p, e, t);

plot(basisobj)

coef = nodes(:,1).^2 + nodes(:,2).^2;

fdobj = fd(coef,basisobj);

plot(fdobj)

data = zeros(nNodes,2);


data(:,1) = (1:nNodes)';
data(:,2) = nodes(:,1).^2 + nodes(:,2).^2;

lambda = 100;
newfdobj = smooth_FEM_fd(data,fdobj,lambda);
plot(newfdobj)

