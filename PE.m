
%% 
clear all,clc,close all
%West Palm
filename = 'Book1.xlsx';
sheet = 'West Palm';
y_true = 'B21:C69';

T=0.25;
% t=0:0.25:12-(0.25*5);
t=0:0.25:12;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('per day')
legend ('\beta_t')
title(sheet)
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('per day')
xlabel('Days')


beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_WP=t(row);
t_r_WP=t(col);
beta_WP1= x(3,:)';
rec_WP1=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_WP,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_WP,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );


figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')


% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% Miami-Ft. Lauderdale
filename = 'Book1.xlsx';
sheet = 'MiamiFtLauderdale';
y_true = 'B21:C69';

T=0.25;
% t=0:0.25:12-(0.25*5);
t=0:0.25:12;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('per day')
legend ('\beta_t')
title(sheet)
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('per day')
xlabel('Days')


beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_MIA=t(row);
t_r_MIA=t(col);
beta_MIA= x(3,:)';
rec_MIA=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_MIA,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_MIA,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );

e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')

% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% FortMyersNaples
filename = 'Book1.xlsx';
sheet = 'FortMyersNaples';
y_true = 'B21:C69';

T=0.25;
% t=0:0.25:12-(0.25*5);
t=0:0.25:12;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);

beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_FTM=t(row);
t_r_FTM=t(col);
beta_FTM= x(3,:)';
rec_FTM=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_FTM,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_FTM,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')

% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% TampaStPetersburg
filename = 'Book1.xlsx';
sheet = 'TampaStPetersburg';
y_true = 'B21:C69';

T=0.25;
% t=0:0.25:12-(0.25*5);
t=0:0.25:12;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('per day')
legend ('\beta_t')
title(sheet)
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('per day')
xlabel('Days')


beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_TPA=t(row);
t_r_TPA=t(col);
beta_TPA= x(3,:)';
rec_TPA=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_TPA,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_TPA,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')
% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% Orlando
filename = 'Book1.xlsx';
sheet = 'Orlando';
y_true = 'B25:C69';

T=0.25;
% t=0:0.25:12-(0.25*5);
t=0:0.25:11;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('per day')
legend ('\beta_t')
title(sheet)
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('per day')
xlabel('Days')


beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_MCO=t(row);
t_r_MCO=t(col);
beta_MCO= x(3,:)';
rec_MCO=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_MCO,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_MCO,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')

% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% Jacksonville
filename = 'Book1.xlsx';
sheet = 'Jacksonville';
y_true = 'B26:C69';

T=0.25;
t=0:0.25:12-(0.25*5);


N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true);
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('per day')
legend ('\beta_t')
title(sheet)
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('per day')
xlabel('Days')


beta= x(3,:);
rec=x(4,:);
b_WP=beta(row);
r_WP=rec(col);
t_b_JAX=t(row);
t_r_JAX=t(col);
beta_JAX= x(3,:)';
rec_JAX=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=0.25;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_WP*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_WP*xt(1,k-1)*xt(2,k-1)-r_WP*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t_b_JAX,beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t_r_JAX,rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true,'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true' S_true' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')

% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')

%% IRMA beta and rec
beta_all=[0.008 0.0111 0.0089 0.01 0.006 0.00787];
t_b=[t_b_WP t_b_MIA t_b_FTM t_b_TPA t_b_MCO t_b_JAX];
rec_all=[0.195056809 0.1841 0.1901 0.1708 0.2214 0.1799];
t_r=[t_r_WP t_r_MIA t_r_FTM t_r_TPA t_r_MCO t_r_JAX];
t=0:0.25:12;t1=0:0.25:11;t2=0:0.25:12-(0.25*5);
figure
subplot(211)
plot(t,beta_WP1,'-k+', t, beta_MIA,'-ro', t, beta_FTM,'-bx', t, beta_TPA,'-g^', t1, beta_MCO,'-md', t2, beta_JAX,'--k')
%% Wilmington
filename = 'Book1.xlsx';
sheet = 'Wilmington';
y_true = 'C2:D445';

T=1/24;

t=0:T:18.4583;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true(1:end-1));
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);

beta= x(3,:);
rec=x(4,:);
b_Wilm=beta(row);
r_Wilm=rec(col);
beta_Wilm= x(3,:)';
rec_Wilm=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=T;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_Wilm*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_Wilm*xt(1,k-1)*xt(2,k-1)-r_Wilm*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t(row),beta(row),'ro')
ylabel('Transmission rate per capita (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t(col),rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true(1:end-1),'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')

legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );


% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}')
e=[t' I_true(1:end-1)' S_true(1:end-1)' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')



%% GNW
filename = 'Book1.xlsx';
sheet = 'GreenvilleNew BernWashington';
y_true = 'C2:D445';

T=1/24;

t=0:T:18.45833333;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true(1:end-1));
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('Transmission rate per capita (\beta)')
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('REcovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );


beta= x(3,:);
rec=x(4,:);
b_GNW=beta(row);
r_GNW=rec(col);
beta_GNW= x(3,:)';
rec_GNW=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=T;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_GNW*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_GNW*xt(1,k-1)*xt(2,k-1)-r_GNW*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t(row),beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t(col),rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true(1:end-1),'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );
e=[t' I_true(1:end-1)' S_true(1:end-1)' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')
% figure(3)
% subplot(2,1,1)
% plot(t,x(1,:),'r',t,S_true,'--bo')
% legend ('S_{est}','S_{true}')
% subplot(2,1,2)
% plot(t,x(2,:),'r',t,I_true,'--bo')
% ylabel('Percentage of Gas Stations')
% legend ('I_{est}','I_{true}'

%% RDF
filename = 'Book1.xlsx';
sheet = 'Raleigh-DurhamFayetteville';
y_true = 'C2:D445';

T=1/24;

t=0:T:18.4583;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true(1:end-1));
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('Transmission rate per capita (\beta)')
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );


beta= x(3,:);
rec=x(4,:);
b_GNW=beta(row);
r_GNW=rec(col);
beta_GNW= x(3,:)';
rec_GNW=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=T;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_GNW*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_GNW*xt(1,k-1)*xt(2,k-1)-r_GNW*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t(row),beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t(col),rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true(1:end-1),'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );

e=[t' I_true(1:end-1)' S_true(1:end-1)' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')
%% NPN
filename = 'Book1.xlsx';
sheet = 'Norfolk-Portsmouth-Newport News';
y_true = 'C4:D445';

T=1/24;

t=0:T:18.4583-2*T;

N = size(t,2); % Number of time steps for filter

% Step 1: Define UT Scaling parameters and weight vectors
L = 4; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
Q = diag([1e1 1e1 1e2 1e2]);
R = diag([1e1 1e1]);


%Step 3: Observation
yt=(xlsread(filename,sheet,y_true))'*100;
S_true= yt(1,:);I_true=yt(2,:);


% Step 4: Initialize state and covariance
x = zeros(L, N); % Initialize size of state estimate for all k
% binit= (yt(2,2)-yt(2,1))/(yt(1,1)*yt(2,1)*T);
% rinit= (yt(3,2)-yt(3,1))/(yt(2,1)*T);
x(:,1) = [yt(:,1);0;0]; % Set initial state estimate
P0 = eye(L,L); % Set initial error covariance

% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0

for k = 2:N

% Step 1: Generate the sigma-points
sP = chol(P); % Calculate square root of error covariance
% chi_p = "chi previous" = chi(k-1)
chi_p = [x(:,k-1) x(:,k-1)*ones(1,L)+sqrt(L+lambda)*sP, ...
x(:,k-1)*ones(1,L)-sqrt(L+lambda)*sP];

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
%chi_m = "chi minus" = chi(k|k-1)
for i=2:2*L+1
chi_m(:,i) =[chi_p(1,i)+((-chi_p(3,i)*chi_p(1,i)*chi_p(2,i)))*T;... //S
           chi_p(2,i)+(chi_p(3,i)*chi_p(1,i)*chi_p(2,i)-chi_p(4,i)*chi_p(2,i))*T;... //I
           chi_p(3,i);... //beta
           chi_p(4,i)]; %r
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q;
for i = 1:2*L+1
P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
psi_m = [(chi_m(1,:));(chi_m(2,:))];

y_m = psi_m*wm; % Calculate mean of predicted output
% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
Pxy = zeros(L,2);
for i = 1:2*L+1
Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate
end
beta_WP=(x(3,:));
r_WP=(x(4,:));

S=zeros(length(beta_WP),length(t),length(r_WP));
I=zeros(length(beta_WP),length(t),length(r_WP));

for m=1:1:length(beta_WP)

    S(m,1,:)=yt(1,1); %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(m,1,:)=yt(2,1);
   
end

for i=1:length(x(3,:))
    for j=1:length(x(4,:))
        for k = 2:N
            S(i,k,j)= S(i,k-1,j)+((-beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)*T));
            I(i,k,j)= I(i,k-1,j)+(beta_WP(i)*S(i,k-1,j)*I(i,k-1,j)-r_WP(j)*I(i,k-1,j))*T; 
        end
    MSE(i,j)=immse(I(i,:,j),I_true(1:end-1));
    end
end

[maxNum, maxIndex] = min(MSE(:));
[row, col] = ind2sub(size(MSE), maxIndex);
figure
subplot(2,1,1)
plot(t,x(3,:))
ylabel('Transmission rate per capita (\beta)')
subplot(2,1,2)
plot(t,x(4,:))
legend ('r_t')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );


beta= x(3,:);
rec=x(4,:);
b_GNW=beta(row);
r_GNW=rec(col);
beta_GNW= x(3,:)';
rec_GNW=x(4,:)';
xt = zeros(2, N); % Initialize size of true state for all k
xt(:,1) = [yt(:,1)]; % Set true initial state
T=T;
cf_WP=1;
for k = 2:N
xt(:,k) = [xt(1,k-1)+((-b_GNW*xt(1,k-1)*xt(2,k-1))*cf_WP*T);... //S
           xt(2,k-1)+(b_GNW*xt(1,k-1)*xt(2,k-1)-r_GNW*xt(2,k-1))*T*cf_WP];
           
end
WestPalm=xt;
t1=0.75:0.25:12.75;

figure
subplot(2,1,1)
plot(t,x(3,:),t(row),beta(row),'ro')
ylabel('Transmission rate Per Captia (\beta)')
subplot(2,1,2)
plot(t,x(4,:),t(col),rec(col),'ro')
ylabel('Recovery rate (\gamma)')
xlabel('Days')
fn = sprintf('Irma_beta_gamma_%s',sheet); 
saveas( gcf, fn, 'png' );

figure
plot(t,xt(2,:),'r',t,I_true(1:end-1),'--b')
xlabel('Days')
ylabel('Percentage of Gas Stations')
% title(sheet)
legend ('I_{UKF}','I_{true}')
fn = sprintf('Irma_%s',sheet); 
saveas( gcf, fn, 'png' );

e=[t' I_true(1:end-1)' S_true(1:end-1)' x(1,:)' x(2,:)', x(3,:)', x(4,:)'];  
xlswrite('Irma Data.xls', e, sheet, 'A2')