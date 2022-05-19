%% This loads our data
[X, y] = load_data_ex2();

%% Normalise and initialize.
[X, mean_vec, std_vec] = normalise_features(X);

%after normalising we add the bias
X = [ones(size(X, 1), 1), X];

%initialise theta
theta = [0.0, 0.0, 0.0];
alpha = 0.1;
iterations = 100;

%% 
t = gradient_descent(X, y, theta, alpha, iterations);
disp 'Press enter to exit!';
pause;

disp(t)
house1Area = (1650-mean_vec(1))/std_vec(1); %modified
house1Bedrooms = (3-mean_vec(2))/std_vec(2); %modified
house1 = [1 house1Area house1Bedrooms]; %modified

house2Area = (3000-mean_vec(1))/std_vec(1);%modified
house2Bedrooms = (4-mean_vec(2))/std_vec(2); %modified
house2 = [1 house2Area house2Bedrooms]; %modified

price1 = house1*t.'; %modified
price2= house2*t.'; %modified

disp(price1);
disp(price2);