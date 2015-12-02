% testSigmoid.m
%
% Author: Wen Shinhua
%
% Modified by Antonio Ulloa on Tue Jul  8 12:57:06 EDT 2003
%
% July 8th, 2003

function testSigmoid

MAX=243;                     % maximum number of units

BETA1=1;                     % BETA controls sigmoid's upper limit
BETA2=1;

SIGMA1=10 ;                  % SIGMA controls the rate of rise of sigmoid: 
SIGMA2=100;                  % higher values of SIGMA slow down sigmoid

THETA1=18;                   % threshold for sigmoid activation
THETA2=0;

r=0:MAX;                     % creates a vector containing [1..243]

match=sigmoid(r,BETA1,SIGMA1, THETA1);

nonMatch=sigmoid(r,BETA2,SIGMA2, THETA2);

plot(r,match, r, nonMatch); 

% ========================= function y=sigmoid(r) ==============
function fx=sigmoid(x,BETA,SIGMA, THETA)   % x may be a vector 

fx=BETA.*pos(x-THETA)./(pos(x-THETA)+SIGMA);

% @@@@@@@@@@@@@@@@@@@@@@@@@ End function y=sigmoid(r) @@@@@@@@@@

