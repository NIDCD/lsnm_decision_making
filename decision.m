% Takes a decision on a delayed match-to-sample experiment
%
% Antonio Ulloa (Based on Shihua Wen's code)
% Brain Imaging & Modeling Section, NIDCD/NIH
%
% Tue Jun 17 13:37:12 EDT 2003
%
% Last updated: Tue Jul  8 14:44:22 EDT 2003

%--------------------------------------------------------------------------
% Intialize parameters
%--------------------------------------------------------------------------
MAX=27;                     % maximum number of active response units

t0=0;                       % initial time for integration
tf=3.67;                    % final time for integration
dt=0.01;                    % time step
tSpan=[t0:dt:tf];           % time span

global R;                   % dynamic vector with # of active units

R1=load('ABCvsABCexfr1.out');
R2=load('ABCvsABCexfr2.out');
R3=load('ABCvsABCexfr3.out');

R=[R1 R2 R3];               % Concatenation of the response matrices

[t,x] = ode45('circuit', [tSpan], [0, 0, 0, 0, 0] );

figure;
plot(t, x(:,1), 'r', t, x(:,2),'b');
legend('match', 'nonmatch');

figure;
plot(t, x(:,4), 'r', t, x(:,5), 'b');
legend('match sigmoid', 'nonmatch sigmoid');

figure;
plot(t, x(:,3), 'k');
legend('active units');
