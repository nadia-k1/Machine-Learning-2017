%%
%load phoneme1
%phoneme1 = load('modelPhoneme1.mat');%model with k=3
 phoneme1 = load('modelPhoneme1B.mat');%model with k=6
phoneme1_mu = phoneme1.mu;
phoneme1_s2 = phoneme1.s2;
phoneme1_p = phoneme1.p;

%%
%load phoneme2
%phoneme2 = load('modelPhoneme2.mat');%model with k=3
 phoneme2 = load('modelPhoneme2B.mat');%model with k=6
phoneme2_mu = phoneme2.mu;
phoneme2_s2 = phoneme2.s2;
phoneme2_p = phoneme2.p;

%%
dataVector1 = test_set1(:,2:end);
dataVector2 = test_set2(:,2:end);
testData = [dataVector1 ; dataVector2];
[m d] = size(testData);
actual_class_k6 = [test_set1(:,1);test_set2(:,1)];
k=6;

%% calculate probability 
ml1 = 0.0;
for i=1:k
    ml1 = ml1+ (phoneme1_p(i)/(2*pi*det(phoneme1_s2(:,:,i)))^(0.5))*exp(-0.5*sum((testData'-repmat(phoneme1_mu(:,i),1,m))'*inv(phoneme1_s2(:,:,i)).*(testData'-repmat(phoneme1_mu(:,i),1,m))',2));
end

ml2 = 0.0;
for i=1:k
    ml2 = ml2+ (phoneme2_p(i)/(2*pi*det(phoneme2_s2(:,:,i)))^(0.5))*exp(-0.5*sum((testData'-repmat(phoneme2_mu(:,i),1,m))'*inv(phoneme2_s2(:,:,i)).*(testData'-repmat(phoneme2_mu(:,i),1,m))',2));
end

probabilities_k6 = [ml1 ml2];
% save('probabilities_k6.mat', 'probabilities_k6');

%%
pred_class_k6 = zeros(length(testData),1);
pred_class_k6(probabilities_k6(:,1)>probabilities_k6(:,2)) = 1;%assign class 1 if probability for phoneme 1 > probability for phoneme 2
pred_class_k6(probabilities_k6(:,1)<probabilities_k6(:,2)) = 2;%assign class 2 if probability for phoneme 1 < probability for phoneme 2
% save('pred_class_k6.mat', 'pred_class_k6');

cfmat1 = pred_class_k6(:,1)==actual_class_k6(:,1); %compare if actual outputs matches with prediction
% save('cfmat1.mat', 'cfmat1');
% correct = (sum(pred_class_k6(:,1)==actual_class_k6(:,1))/length(testData))*100;
error_k6 = sum(cfmat1(:,1)==0)/length(testData); %calculate miss classified predictions
display (error_k6);
% save('error_k6.mat', 'error_k6');