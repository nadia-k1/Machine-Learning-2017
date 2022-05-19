%% create grid of points
f1_total = [X1(:,1) ; X2(:,1)]; %combine together f1 for phoneme 1 and 2
f2_total = [X1(:,2) ; X2(:,2)]; %combine together f2 for phoneme 1 and 2
x1 = min(f1_total):max(f1_total); %create vector with range of minimum and maximum of f1 for phoneme 1 and 2
y1 = min(f2_total):max(f2_total); %create vector with range of minimum and maximum of f2 for phoneme 1 and 2
[F1,F2] = meshgrid(x1, y1); 
points = [F1(:) F2(:)]; %create vector of points with all possible combinations
[z l] = size(points)

M = [];
k=3;

%% calculate probabilities
ml1 = 0.0;
        for a=1:k
            ml1 = ml1+ (phoneme1_p(a)/(2*pi*det(phoneme1_s2(:,:,a))^(0.5)))*exp(-0.5*sum((points'-repmat(phoneme1_mu(:,a),1,z))'*inv(phoneme1_s2(:,:,a)).*(points'-repmat(phoneme1_mu(:,a),1,z))',2));
        end

ml2 = 0.0;
        for b=1:k
            ml2 = ml2+ (phoneme2_p(b)/(2*pi*det(phoneme2_s2(:,:,b))^(0.5)))*exp(-0.5*sum((points'-repmat(phoneme2_mu(:,b),1,z))'*inv(phoneme2_s2(:,:,b)).*(points'-repmat(phoneme2_mu(:,b),1,z))',2));
        end

probabilities = [ml1 ml2];

%% create classification vector
classification(probabilities(:,1)>=probabilities(:,2)) = 1;
classification(probabilities(:,1)<probabilities(:,2)) = 2;

%% create matrix with the predicted classes
M = reshape(classification,length(y1),length(x1));

%% plot matrix
imagesc(M);


