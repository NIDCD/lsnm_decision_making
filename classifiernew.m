% The Competitive learning classifer
% Author: Shihua Wen
% DATE: June-16-2003
% UPDATED: Sep-22-2003
%disp('I am here');
% ======================== Main function Decision = classifier() ==========================
function Decision = classifier
global Cha ActiveUnits DataSeq Totrials TimeSteps;

% Load or not load new file
choi='N';   % input('Need loading New data file? [Y/N]: ','s');
switch choi
    case {'Y','y'}
       [ActiveUnits,Cha,Totrials,DataSeq,TimeSteps,...
        NumofUnits,Threshold]=DATAprepare; 
       disp('New Data loaded ! ');
   case {'N','n'} 
       [ActiveUnits,Cha,Totrials,DataSeq,TimeSteps,...
        NumofUnits,Threshold]=loadata;
       disp('Old Data loaded ! ');
    otherwise
       disp('Input wrong, program stopped!!!');    
end

%speedecide;

disp('Which sound you want to compare with sound ABC?');
Aseq=input(' Type 1~28 :','s');   
rtrialseq=str2num(Aseq);
info4disp=strcat('The sound you select is: ',Cha{1,rtrialseq});
disp(info4disp);

% Parameters for the competitive system equation
options = odeset('RelTol',1e-3,'AbsTol',[1e-3 1e-3 1e-3]);
Xinit=[0 0 1];  % Initial start
tspan=[0:0.02:7.36];   % Duration time of the input signal
% Parameters of Competitive learning model and its sigmoid function 
% Fixed parameters not varied by decision 
rA=[1 1 0.3]';     % Decay rate 
rB=[1 1 1]';     % Upper activity bound   
rC=[0 0 0]';     % lower activity bound

% various parameters varied by decision [% 1:Match , 2:Non-mathch, 3: stay]
rTheta=[0 0 0]';    % Summed activity threshold 
rBeta=[2 2 1]';     % Gain factor
rSigma=[0.5 10 10]';    % Sigmoidal threshold
rz=0.02;            % 20ms in each time step


% Call core program to solve the competitive system
[Ttotal,Xtotal]=ode45(@comp,tspan,Xinit,options,rA,rB,rC,rTheta,rBeta,rSigma,rz,rtrialseq); 

%Plotting the columns of the returned array Y versus T shows the solution 
figname=strcat(Cha{1,rtrialseq},' :  blue is X(1)match; Red is X(2)notmatch; Black is X(3)stay');
figure('name',figname,'NumberTitle','off');
plot(Ttotal,Xtotal(:,1),'-',Ttotal,Xtotal(:,2),'-r',Ttotal,Xtotal(:,3),'-k');
%xlim([0 ceil(timeinterval)]);
ylim([0 1]);


disp('Program finished');


% @@@@@@@@@@@@@@@@@@@@@@@@@ End Main function Decision=classifier()  @@@@@@@@@@@@@@@@@@@@@@@@





% ========================= function dX=comp(t,X) =========================
% ===== The Competitive network =====
function dX=comp(t,X,A,B,C,Theta,Beta,Sigma,z,trialseq)
global TimeSteps;

if fix(t/z)+1<=TimeSteps
   sts=fix(t/z)+1;
else   
   sts=TimeSteps;
end   

see=[t,sts]

r=DATAinputC(trialseq,sts); % Activity of the response units from NIH;
%rneuron(1)=r;
%rneuron(2)=r*2/6;
%rneuron(3)=0;

dX=zeros(3,1);    % initial value of Xs.  X(1) is mathch, X(2) is non-mathch, X(3) is stay  

% follow pepe's draft
dX(1)=-A(1)*X(1)+(B(1)-X(1))*(sigmoid(r,Beta(1),Sigma(1)))-(C(1)+X(1))*(X(2)+X(3));
dX(2)=-A(2)*X(2)+(B(2)-X(2))*(1-sigmoid(r,Beta(2),Sigma(2)))-(C(2)+X(2))*(X(1)+X(3));
dX(3)=-A(3)*X(3)+(B(3)-X(3))*(sigmoid(0,Beta(3),Sigma(3)))-(C(3)+X(3))*(X(1)+X(2));

% my combined one
%dX(1)=-A(1)*X(1)+(B(1)-X(1))*(sigmoid(r,Beta(1),Sigma(1)));
%dX(2)=-A(2)*X(2)+(B(2)-X(2))*(1-sigmoid(r,Beta(2),Sigma(2)));
%dX(3)=-A(3)*X(3)+(B(3)-X(3))*(sigmoid(0,Beta(3),Sigma(3)));

% @@@@@@@@@@@@@@@@@@@@@@@@@ End function dX=comp(t,x) @@@@@@@@@@@@@@@@@@@@@@@@





% ========================= function y=sigmoid(r) =========================
function y=sigmoid(z,vBeta,vSigma)   % z may be a vector 

%y=vBeta.*(z-vTheta)./((z-vTheta)+vSigma);
y=vBeta.*z./(z+vSigma);
%y=vBeta./exp(-(z-vTheta)./vSigma);
%y=exp(-(z-vSigma)./vBeta);

% @@@@@@@@@@@@@@@@@@@@@@@@@ End function y=sigmoid(r) @@@@@@@@@@@@@@@@@@@@@@@@@




% ============================== function Datainput for more (density) =====================
function TotSmore=DATAinputD(vtrialseq,vsts)    % TotSmore is only one number meaning how many 
                                               % ActiveUnits in the current time step are more
                                               % than the number in the total previous steps. 
% recieve data from the signal
global Cha DataSeq Totrials NumofUnits Threshold;

Singalin=zeros(1,3);
Smore=zeros(1,3);
for UnitNo=1:3 
   Worksingal=DataSeq{vtrialseq,UnitNo}(1:vsts,:);
   for internalUnit=1:NumofUnits
         if length(find(Worksingal(:,internalUnit)>=Threshold))>0
             Actornot(internalUnit)=1;
         else
             Actornot(internalUnit)=0;
         end   
   end 
   if sum(Actornot)>Singalin(UnitNo)
       Smore(UnitNo)=sum(Actornot)-Singalin(UnitNo); 
   else    
       Smore(UnitNo)=0;
   end   
   Singalin(UnitNo)=sum(Actornot);
   % pause(0.02);
end
TotSmore=sum(Smore);

% @@@@@@@@@@@@@@@@@@@@@@@@@ End function Datainput for more (density) @@@@@@@@@@@@@@@@@@@@@@@@@





% ============================== function Datainput Cumulative =====================
function TotSingal=DATAinputC(vtrialseq,vsts)    % TotSmore is only one number meaning how many 
                                               % ActiveUnits in the current time step are more
                                               % than the number in the total previous steps. 
% recieve data from the signal
global Cha DataSeq Totrials NumofUnits Threshold;

Singalin=zeros(1,3);
for UnitNo=1:3 
   Worksingal=DataSeq{vtrialseq,UnitNo}(1:vsts,:);
   for internalUnit=1:NumofUnits
         if length(find(Worksingal(:,internalUnit)>=Threshold))>0
             Actornot(internalUnit)=1;
         else
             Actornot(internalUnit)=0;
         end   
   end 
   Singalin(UnitNo)=sum(Actornot);
   % pause(0.02);
end
TotSingal=sum(Singalin);

% @@@@@@@@@@@@@@@@@@@@@@@@@ End function Datainput Cumulative @@@@@@@@@@@@@@@@@@@@@@@@@





% ============================== function Dataprepare =====================
function [ActiveUnits,Cha,Totrials,DataSeq,totlocs,...
          TimeSteps,NumofUnits,Threshold]=DATAprepare
% Counting Active Units in each response units in each trials 

global ActiveUnits Cha Totrials DataSeq totlocs...
       TimeSteps NumofUnits Threshold;
      
load MatchornotDataset;
ActiveUnits=zeros(28,3);   % Counting active units which are greater than the threshold 

TimeSteps=368;
NumofUnits=81;
Threshold=0.6;
LocUnit=1;
% Index table for reference
% Cha={'AAA' 'AAB' 'AAC'       'ABA' 'ABB' 'ABC'      'ACA' 'ACB' 'ACC' ...
%     'BAA' 'BAB' 'BAC'       'BBA' 'BBB' 'BBC'      'BCA' 'BCB' 'BCC' ...
%     'CAA' 'CAB' 'CAC'       'CBA' 'CBB' 'CBC'      'CCA' 'CCB' 'CCC' ...
%     'DBF'};

for trial=1:Totrials   % Totrials=28
   for responseneuron=1:3
      WorkingSet=DataSeq{trial,responseneuron};   
      Actornot=zeros(NumofUnits,1);
      for UnitNo=1:NumofUnits
         Locs=find(WorkingSet(:,UnitNo)>=Threshold);
         LL=length(Locs);
         if LL>0
             Actornot(UnitNo)=1;
             totlocs(LocUnit,:)=[trial,responseneuron,UnitNo,Locs(1),Locs(LL)];  % Units, first time>threshold, last time>threshold      
             LocUnit=LocUnit+1;
             %disp('I have been here.');
         else
             Actornot(UnitNo)=0;
         end    
         %plot(ABCvsAACexfr1(:,i));
      end 
      ActiveUnits(trial,responseneuron)=sum(Actornot); 
      
   end
   funsavar=sprintf('%s=ActiveUnits(%d,:);',strcat('ABCvs',Cha{trial}),trial);
   eval(funsavar); 
  
end

save('Counting','ActiveUnits','Cha','Totrials','DataSeq','totlocs',...
                'TimeSteps','NumofUnits','Threshold');    
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ End function Dataprepare @@@@@@@@@@@@@@@@@@@@@@@@



% ================== Function of loading old data =================
function [ActiveUnits,Cha,Totrials,DataSeq,TimeSteps,...
        NumofUnits,Threshold]=loadata
global ActiveUnits Cha Totrials TimeSteps;
load Counting;
% @@@@@@@@@@@@@@@@@@@ End function of loading data @@@@@@@@@@@@@@@@@




% Original put the
%********************************* function speedecide *********************************
function speedecide
global Cha ActiveUnits Totrials;

contin='C';
while contin=='C' | contin=='c'    
disp('Which sound you want to compare with sound ABC?');
seq=input(' Type 1~28 :','s');
info4disp=strcat('The sound you select is: ',Cha{1,str2num(seq)});
disp(info4disp);
match1=0;
match2=0;
match3=0;
if ActiveUnits(str2num(seq),1)>=5
    match1=1;
    if ActiveUnits(str2num(seq),2)>=5
       match2=1;
       if ActiveUnits(str2num(seq),3)>=5    
           match3=1;
       else
           match3=0;
       end    
    else  
       match2=0; 
    end   
else   
    match1=0;  
end

if (match1*match2*match3==1)
   decision=1;
   disp('Result is:   Match');
else
   decision=2;
   disp('Result is:   Not match');
end   

disp('Continue testing?');
contin=input('Press ''C'' to continue, or press any ''Ctrl+C'' to stop: ','s');
end

% ****************** End function speedecide ******************************************************