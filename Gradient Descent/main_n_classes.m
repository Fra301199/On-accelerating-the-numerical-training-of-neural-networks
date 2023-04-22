%% PRE-PROCESSING

clear all; close all; clc;

% synthetic dataset generated
load X.mat
load y.mat

% normalize dataset
m = mean(X);
s = max(X);
for i = 1:2
    X(:,i) = (X(:,i) + m(i))/s(i);   % normalize every dimension
end

% parameters of the dataset
n = size(X, 1);
d = size(X, 2);
n_classes = 4;

figure
scatter(X(:, 1), X(:, 2), [], y);

% transform y to categorical (inefficient but needed)
y_cat = zeros(4, n);
for k = 1:n
    y_cat(y(k)+1,k) = 1;
end

% decide how to split the dataset
split_ratio = 0.3;
k = floor(split_ratio*n);

% split dataset into training and test
X_train = X(1:k,:);
X_test = X(k+1:end,:);
y_train = y(1:k);
y_test = y(k+1:end);
y_train_cat = y_cat(:,1:k);
y_test_cat = y_cat(:,k+1:end);

% define activation functions
sigma = @(t) 1./(1+exp(-t));
sigmaprime = @(t) sigma(t).*(1-sigma(t));

% define hyperparameters of the network ANY NETWORK SHAPE (this includes 
% input & output layers too)
shape = [d, 3, 3, n_classes];

niter = 3e5;    % max iterations
eta = 0.05;     % learning rate

%% TRAINING
[costHistory, W, b] = GradientDescent( ...
        X_train', y_train_cat, niter, sigma, sigmaprime, eta, shape);

% save the network
save NNparams.mat W b sigma sigmaprime

% plot training behaviour
figure
plot(linspace(0,niter,niter)', costHistory, '-')
fprintf('Cost Function: %f\n', costHistory(end));

filename = 'training_cost_function.png'; 
saveas(gcf, filename, 'png');
%% TEST
load NNparams.mat W b

% predict the classes of the test points using the network
y_pred = Predict(W, b, sigma, X_test');

% plot predictions
figure
scatter(X_train(:, 1), X_train(:, 2), [], y_train, 'o');
hold on;
scatter(X_test(:, 1), X_test(:, 2), [], y_pred, '*');
saveas(gcf, '4_classes.png', 'png');

% confusion matrix
figure
confusionchart(y_test,y_pred)
saveas(gcf, 'ConfusionMatrix', 'png');

%% 
load NNparams.mat W b
% generate new random data in the square [0,1]^2
n = 1000;
testDatax = -1 + 2*rand(n, 2);
testDatay = Predict(W, b, sigma, testDatax');

figure
scatter(X_train(:,1), X_train(:,2), [], y_train, 'o');
hold on;
scatter(testDatax(:,1), testDatax(:,2), [], testDatay, 'filled');
saveas(gcf, 'SpaceExploration', 'png');