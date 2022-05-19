% Simple script to do EM for a mixture of Gaussians.
% -------------------------------------------------
%  based on code from  Rasmussen and Ghahramani
% (http://www.gatsby.ucl.ac.uk/~zoubin/course02/)



% load('PB12.mat');
%% phoneme1
% scrambled_indexes = randperm(size(X1,1));
% training_set1=X1(scrambled_indexes(1:106),:);
% test_set1=[ ones(46,1) X1(scrambled_indexes(107:end),:)];
% save('training_set1.mat','training_set1');
% save('test_set1.mat','test_set1');

%% phoneme2
% scrambled_indexes1 = randperm(size(X2,1));
% training_set2=X2(scrambled_indexes1(1:106),:);
% test_set2=[ 2*ones(46,1) X2(scrambled_indexes(107:end),:)];
% save('training_set2.mat', 'training_set2');
% save('test_set2.mat', 'test_set2');

%% % Initialise parameters
tr1 = load('training_set1.mat');
x = [tr1.training_set1 tr1.training_set1(:,1)+tr1.training_set1(:,2)];
% x = load('training_set2');
[n D] = size(x);                    % number of observations (n) and dimension (D)
k = 3;                              % number of components
p = ones(1,k)/k;                    % mixing proportions
mu = x(ceil(n.*rand(1,k)),:)';      % means picked rqandomly from data
s2 = zeros(D,D,k);                  % covariance matrices
niter=100;                          % number of iterations

% initialize covariances 

   disp(cov(x)./3);
for i=1:k
  s2(:,:,i) = cov(x)./k;      % initially set to fraction of data covariance
end

set(gcf,'Renderer','zbuffer');

clear Z;
try

  % run EM for niter iterations
  
  for t=1:niter,
    fprintf('t=%d\r',t);
    % Do the E-step:
  
    for i=1:k
      Z(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));
    end
    Z = Z./repmat(sum(Z,2),1,k);
    
    % Do the M-step:
    
    for i=1:k
      mu(:,i) = (x'*Z(:,i))./sum(Z(:,i));
      
      % We will fit Gaussians with diagonal covariances:
      
%       s2(:,:,i) = diag((x'-repmat(mu(:,i),1,n)).^2*Z(:,i)./sum(Z(:,i))); 
      
      % To fit general Gaussians use the line:
      diag_mat = eye(D)*0.01;
      
      s2(:,:,i) = (x'-repmat(mu(:,i),1,n))*(repmat(Z(:,i),1,D).*(x'-repmat(mu(:,i),1,n))')./sum(Z(:,i))+diag_mat;
      
      p(i) = mean(Z(:,i));
    end    
    
    clf
    hold on
    plot(x(:,1),x(:,2),'.');
    for i=1:k
      plot_gaussian(2*s2(:,:,i),mu(:,i),i,11);
    end
    drawnow;
  end
%   save('modelPhoneme1B.mat', 'mu', 's2', 'p');
%   save('modelPhoneme2B.mat', 'mu', 's2', 'p');
catch
  disp('Numerical Error in Loop - Possibly Singular Matrix');
end;