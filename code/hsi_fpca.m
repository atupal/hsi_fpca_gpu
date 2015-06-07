% 基于函数型主成分分析的高光谱图像分类

close all;
clear;

% 导入数据
load Indian.mat;
% load KSC.mat;
% load Botswana.mat;
% load Pavia.mat;

[train_data,test_data] = scaleForSVM(train_data,test_data,0,1);
[numtrain,dimtrain] = size(train_data);   %训练样本大小和维数
[numtest,dimtest] = size(test_data);      %测试样本大小和维数
sample = [train_data', test_data'];
num = numtrain + numtest;

r = 10^(-3);
rangewavelength = r * [1, dimtrain]';
norder = 4;                % B样条阶数
wavelength = r * (1:dimtrain)';   % 波长
nbasis = norder + length(wavelength) - 2;
Lfd = 2;

% 生成B样条基函数
basisobj = create_bspline_basis(rangewavelength, nbasis, norder, wavelength);

fdnames{1} = 'Wavelength';
fdnames{2} = 'Substance';
fdnames{3} = 'Radiance';

% 选取lambda
lnlam = -8:1:0;
gcvsave = zeros(length(lnlam),1);
for i=1:length(lnlam)
    fdParobj = fdPar(basisobj, Lfd, 10^lnlam(i));
    [~, ~, gcv] = smooth_basis(wavelength, sample, fdParobj, [], fdnames);
    gcvsave(i) = sum(gcv);
end
[~, k] = max(-gcvsave);
lambda = 10^lnlam(k);  % 最优lambda

% 粗糙惩罚法
fdParobj = fdPar(basisobj, Lfd, lambda);
[fdobj, df, gcv] = smooth_basis(wavelength, sample, fdParobj, [], fdnames);

for j = 1:1:40  %函数主成分数目
    
    pcaf = pca_fd(fdobj, j);   %函数型PCA
    
    train_final = pcaf.harmscr(1:numtrain, :);  %训练样本的函数型特征
    test_final = pcaf.harmscr(numtrain + 1:end, :);   %测试样本的函数型特征
    
    %svm分类
    [bestCVaccuracy, bestc, bestg] = SVMcgForClass(train_data_labels,train_final,-4,4,-4,4,10,0.5,0.5);
    cmd = ['-c ', num2str(bestc), ' -g ', num2str(bestg), ' -t ', num2str(2)];
    model = svmtrain(train_data_labels, train_final, cmd);
    type = 1;
    CR = ClassResult(train_data_labels, train_final, model, type);
    type = 2;
    CR2 = ClassResult(test_data_labels, test_final, model, type);
    acc(j) = CR2.accuracy(1);  %精度
    bestgg(j) = bestg;  %最优核参数值
end

