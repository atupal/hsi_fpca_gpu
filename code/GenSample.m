% ����Indian Pines��University of Paviaͼ���ѵ�������Ͳ�������������Ӧ�������

clear

%% ����Indian Pinesͼ��
load Indian_pines_corrected
load Indian_pines_gt

% ����University of Paviaͼ��
% load PaviaU
% load PaviaU_gt

Sample = reshape(indian_pines_corrected,145*145,200);
SampleLabel = reshape(indian_pines_gt,145*145,1);

% Sample = reshape(paviaU, 610*340, 103);
% SampleLabel = double(reshape(paviaU_gt, 610*340, 1));

%% ѡȡ����
ratio = 0.1;   % ѵ����������

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

% �洢Indian Pinesͼ���ѵ�������Ͳ�������������Ӧ�������
save Indian train_data test_data train_data_labels test_data_labels
% �洢University of Paviaͼ���ѵ�������Ͳ�������������Ӧ�������
% save Pavia train_data test_data train_data_labels test_data_labels
