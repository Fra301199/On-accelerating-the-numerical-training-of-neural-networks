clear all; close all; clc;

% visualize data (in this easy example)
x1 = [0.1,0.3,0.1,0.6,0.4];
y1 = [0.1,0.4,0.5,0.9,0.2];

x2 = [0.6,0.5,0.9,0.4,0.7];
y2 = [0.3,0.6,0.2,0.4,0.6];

scatter(x1, y1, 'blue')
hold on;
scatter(x2, y2, 'red')

x = [x1, x2; y1, y2];
y = [ones(1,5) zeros(1,5); zeros(1,5) ones(1,5)];


% define activation functions
sigma = @(t) 1./(1+exp(-t));
sigmaprime = @(t) sigma(t).*(1-sigma(t));

shape = [2, 3, 3, 2];   % ANY NETWORK SHAPE (this includes
%                           input & output layers too)
niter = 3e5;
eta = 0.05;

[costHistory, W, b] = GradientDescent( ...
        x, y, niter, sigma, sigmaprime, eta, shape);

save NNparams.mat W b sigma sigmaprime

% plot results
figure
plot(linspace(0,niter,niter)', costHistory, '-')
fprintf('Cost Function: %f\n', costHistory(end));

%% prediction
load NNparams.mat
% generate new random data in the square [0,1]^2
n = 300;
testDatax = rand(n, 2);
testDatay = Predict(W, b, sigma, testDatax');

x = x';
y = int32([zeros(5, 1); ones(5, 1)]');
figure
scatter(x(:,1), x(:,2), [], y, 'o');
hold on;
scatter(testDatax(:,1), testDatax(:,2), [], testDatay, 'filled');