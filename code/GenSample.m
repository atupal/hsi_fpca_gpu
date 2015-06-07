% 生成Indian Pines或University of Pavia图像的训练样本和测试样本及其相应样本标记

clear

%% 导入Indian Pines图像
load Indian_pines_corrected
load Indian_pines_gt

% 导入University of Pavia图像
% load PaviaU
% load PaviaU_gt

Sample = reshape(indian_pines_corrected,145*145,200);
SampleLabel = reshape(indian_pines_gt,145*145,1);

% Sample = reshape(paviaU, 610*340, 103);
% SampleLabel = double(reshape(paviaU_gt, 610*340, 1));

%% 选取样本
ratio = 0.1;   % 训练样本比例

%%
train_data=[];test_data=[];
train_data_labels=[];test_data_labels=[];

for i = 1:max(SampleLabel)
    ind = find(SampleLabel==i);
    numberofdata(i) = length(ind);
    if(numberofdata(i) ~= 0)
        No = randperm(numberofdata(i));
        data = Sample(ind,:);
        Label = SampleLabel(ind,:);
        
        data1 = data(No(1:ceil(numberofdata(i)*ratio)),:);
        label1 = Label(No(1:ceil(numberofdata(i)*ratio)),:);
        data2 = data(No(ceil(numberofdata(i)*ratio)+1:numberofdata(i)),:);
        label2 = Label(No(ceil(numberofdata(i)*ratio)+1:numberofdata(i)),:);
        
        train_data = [train_data;data1];
        test_data = [test_data;data2];
        train_data_labels = [train_data_labels;label1];
        test_data_labels = [test_data_labels;label2];
    end
end

% 存储Indian Pines图像的训练样本和测试样本及其相应样本标记
save Indian train_data test_data train_data_labels test_data_labels
% 存储University of Pavia图像的训练样本和测试样本及其相应样本标记
% save Pavia train_data test_data train_data_labels test_data_labels
