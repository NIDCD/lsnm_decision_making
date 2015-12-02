% circuit.m
% Equations for decision model
% Antonio Ulloa (Based on Shihua Wen's code)
% Brain Imaging & Modeling Section
% 
% Last updated: Tue Jun 17 13:42:25 EDT 2003

function dx = circuit(t, x)

%--------------------------------------------------------------------------
% Intialize parameters
%--------------------------------------------------------------------------
A=1;           % decay rate of competitive neurons
B=1;           % Upper activity bound of competitive neurons
C=1;           % Lower activity bound of competitive neurons

global R;      % matrix contaning all response units

BETA_m=20;     % BETA controls sigmoid's upper limit     
BETA_n=1;      % BETA controls sigmoid's upper limit
THETA_m=0.18;  % threshold for sigmoid activation
THETA_n=0;     % threshold for sigmoid activation
SIGMA_m=5;     % SIGMA controls the rate of rise of sigmoid:
SIGMA_n=100;   % higher values of SIGMA slow down sigmoid

% The following are sigmoid functions for match (m) and nonmatch (n)
f_m=pos(BETA_m.*(x(3)-THETA_m)./((x(3)-THETA_m)+SIGMA_m));
f_n=pos(BETA_n.*((1-x(3))-THETA_n)./(((1-x(3))-THETA_n)+SIGMA_n));

%--------------------------------------------------------------------------
% Intialize the vector used for derivatives
%--------------------------------------------------------------------------
dx = zeros(5, 1);

%--------------------------------------------------------------------------
% m, match neuron 
%--------------------------------------------------------------------------
dx(1) = -A.*x(1) + (B-x(1)).*f_m - (C + x(1)).*pos(x(2)); 

%--------------------------------------------------------------------------
% n, non-match neuron 
%--------------------------------------------------------------------------
dx(2) = -A.*x(2) + (B-x(2)).*f_n - (C + x(2)).*pos(x(1));

%--------------------------------------------------------------------------
% integrator neuron that counts number of active units in response matrix 
%--------------------------------------------------------------------------
dx(3) = sum(pos(R(round(t.*100+1),:)-0.5));

%--------------------------------------------------------------------------
% This only follows a sigmoid (for plotting purposes)
%--------------------------------------------------------------------------
dx(4) = -x(4) + f_m;

%--------------------------------------------------------------------------
% This only follows a sigmoid (for plotting purposes)
%--------------------------------------------------------------------------
dx(5) = -x(5) + f_n;
