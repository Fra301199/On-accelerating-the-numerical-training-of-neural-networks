% Parareal: first implementation attempt
clear all
close all
clc

% trivial use of parfor

w = [];
tic
for i = 1:30000000
    w(i) = i;
end
toc
% Elapsed time is higher

v = [];
tic
parfor i = 1:30000000
    v(i) = i;
    %v(i+1) = 5 + v(i); not allowed
end
toc
% Elapsed time is lower

% it works on my pc