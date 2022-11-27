clc
clear;
%% default plot settings
set(groot,'defaultLineLineWidth',2.5) 
set(0,'DefaultaxesLineWidth', 1.5) 
set(0,'DefaultaxesFontSize', 14) 
set(0,'DefaultaxesFontWeight', 'bold') 
set(0,'DefaultTextFontSize', 14) 
set(0,'DefaultaxesFontName', 'latex') 
set(0,'DefaultlegendFontName', 'latex')
set(0,'defaultAxesXGrid','off') 
set(0,'defaultAxesYGrid','off') 
set(0,'DefaultTextInterpreter','latex') 
%% read data
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = [10, Inf];
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["TimeSeconds", "Position", "Flag", "Vx_0", "Vy_0", "Vz_0", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"];
opts.SelectedVariableNames = ["TimeSeconds", "Vx_0", "Vy_0", "Vz_0"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];


% Import the data
file = "E:\experiments\experimental data\velocity_data\downstream\y_12x_-3\-3_12_73.Vf";
[filepath,name,ext] = fileparts(file);
A = readmatrix(file,opts);
C = strsplit(name,'_');
x=str2double(C(:,1)); y=str2double(C(:,2)); %z=10;
A_out=fillmissing(A,"linear");
%% Turbulence
%horizontal
u=A_out(:,2);  v=A_out(:,3);  w=A_out(:,4);
umean = mean(u);    vmean = mean(v);    wmean = mean(w);
u_Tseries = u - umean;  v_Tseries = v - vmean; w_Tseries = w - wmean;

sample_variance_u = (length(u_Tseries)-1)*var(u_Tseries)/length(u_Tseries);
sample_variance_v = (length(v_Tseries)-1)*var(v_Tseries)/length(v_Tseries); % estimation of population
sample_variance_w = (length(w_Tseries)-1)*var(w_Tseries)/length(w_Tseries); 

urms=sqrt(sample_variance_u); vrms=sqrt(sample_variance_v); wrms=sqrt(sample_variance_w);

% cross moments (reynolds stresses)
uwmean = mean(u_Tseries.*w_Tseries);
RSS = abs(- 1000*uwmean);  % Reynolds shear stress
uwcorr = uwmean/(sample_variance_w * sample_variance_u)^0.5;   % correlation coef

uvmean = mean(u_Tseries.*v_Tseries); vwmean = mean(v_Tseries.*w_Tseries);   
TKE = 0.5*( sample_variance_u + sample_variance_v + sample_variance_w );
%% higher order moments
%TKE flux
fku=0.5*(mean(u_Tseries.^3)+mean(u_Tseries .* v_Tseries.^2)+mean(u_Tseries .* w_Tseries.^2));
fkv=0.5*(mean(v_Tseries .* u_Tseries.^2)+mean(v_Tseries.^3)+mean(v_Tseries .* w_Tseries.^2));
fkw=0.5*(mean(u_Tseries.^2.*w_Tseries)+mean((w_Tseries .* v_Tseries.^2))+mean(w_Tseries.^3));

% 3rd order
m30 = mean(u_Tseries.^3);
m21 = mean(u_Tseries.^2 .* w_Tseries);
m12 = mean(u_Tseries .* w_Tseries.^2);
m03 = mean(w_Tseries.^3);

skewU = (length(u_Tseries)/(length(u_Tseries)-1))*mean(u_Tseries.^3)/sample_variance_u^(3/2);
skewV = (length(v_Tseries)/(length(v_Tseries)-1))*mean(v_Tseries.^3)/sample_variance_v^(3/2);
skewW = (length(w_Tseries)/(length(w_Tseries)-1))*mean(w_Tseries.^3)/sample_variance_w^(3/2);

% 4th order
m40 = mean(u_Tseries.^4);
m31 = mean(u_Tseries.^3 .* w_Tseries);
m22 = mean(u_Tseries.^2 .* w_Tseries.^2);
m13 = mean(u_Tseries .* w_Tseries.^3 );
m04 = mean(w_Tseries.^4);

kurtU = (length(u_Tseries)/(length(u_Tseries)-1))*mean(u_Tseries.^4)/sample_variance_u^2 - 3;
kurtV = (length(v_Tseries)/(length(v_Tseries)-1))*mean(v_Tseries.^4)/sample_variance_v^2 - 3;
kurtW = (length(w_Tseries)/(length(w_Tseries)-1))*mean(w_Tseries.^4)/sample_variance_w^2 - 3;

%% Quadrants
vel=[u_Tseries w_Tseries];
H=0;
threshold=H*abs(uwmean);
ind=abs(u_Tseries.*w_Tseries)<threshold;
vel(ind,:)=NaN;
q1=find(vel(:,1)>0&vel(:,2)>0); q2=find(vel(:,1)<0&vel(:,2)>0);
q3=find(vel(:,1)<0&vel(:,2)<0); q4=find(vel(:,1)>0&vel(:,2)<0);
q=length(q1)+length(q2)+length(q3)+length(q4);
RS1=length(q1)/q; RS2=length(q2)/q; RS3=length(q3)/q; RS4=length(q4)/q;
%successive events time gap and muliply with 0.04=mean bursting period
e_s=RS2/RS4;
Te=(mean(diff(q2))*0.04);
Ts=(mean(diff(q4))*0.04);
% figure(1)
% plot(u_Tseries, w_Tseries, 'r.')
% hold on
% plot(u_Tseries,abs(uwmean)./u_Tseries)
%% Anisotropy
RST=[mean(u_Tseries.^2)/(2*TKE) - 1/3, mean(u_Tseries.*v_Tseries)/(2*TKE), mean(u_Tseries.*w_Tseries)/(2*TKE) - 1/3
    mean(u_Tseries.*v_Tseries)/(2*TKE), mean(v_Tseries.^2)/(2*TKE) - 1/3, mean(u_Tseries.*w_Tseries)/(2*TKE)
    mean(u_Tseries.*w_Tseries)/(2*TKE), mean(w_Tseries.*v_Tseries)/(2*TKE), mean(w_Tseries.^2)/(2*TKE) - 1/3];
eigen=eig(RST);
I=eigen(1)+eigen(2)+eigen(3);    
II=eigen(1)*eigen(2)+eigen(1)*eigen(3)+eigen(3)*eigen(2) ;     
III=eigen(1)*eigen(2)*eigen(3);
IF=1+27*III+9*II;
Ad=(III/2)/(-II/3)^1.5;
%% output
k=[x y umean vmean wmean urms vrms wrms RSS uwcorr TKE fku fkv fkw m30 m21 m12 m03 skewU skewV skewW m40 m31 m22 m13 m04 kurtU kurtV kurtW RS1 RS2 RS3 RS4 Te Ts IF Ad];
fileID = fopen("tapu1.dat","a+");
fprintf(fileID,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n',k);
%% contours
% N = 500;                           % Number Of Points Desired
% xv = linspace(min(x), max(x), N);
% yv = linspace(min(y), max(y), N);
% [X,Y] = ndgrid(xv, yv);
% Z = griddata(x, y, z, X, Y);
% figure
% contourf(X, Y, Z); shading interp;
% hold on
% rectangle('Position',[1 2 5 6])
% hold on
% rectangle('Position',[1 2 5 6])
% axis('equal')

% %% 2nd time copied
% clc
% clear;
% %% default plot settings
% set(groot,'defaultLineLineWidth',2.5) 
% set(0,'DefaultaxesLineWidth', 1.5) 
% set(0,'DefaultaxesFontSize', 14) 
% set(0,'DefaultaxesFontWeight', 'bold') 
% set(0,'DefaultTextFontSize', 14) 
% set(0,'DefaultaxesFontName', 'latex') 
% set(0,'DefaultlegendFontName', 'latex')
% set(0,'defaultAxesXGrid','off') 
% set(0,'defaultAxesYGrid','off') 
% set(0,'DefaultTextInterpreter','latex') 
% %% fg
% P = 'F:\velocity_confluence\0.2Qr\bed';
% S = dir(fullfile(P,'**','*.dat'));
% for k = 1:numel(S)
%    file = fullfile(S(k).folder,S(k).name);
%    [filepath,name,ext] = fileparts(file);
%     C = strsplit(name,'_');
% end
% %% read data
% function read
% file = "F:\velocity_confluence\0.2Qr\bed\4.4_0.95.dat";
% [filepath,name,ext] = fileparts(file);
% A = load(file);
% C = strsplit(name,'_');
% x=str2double(C(:,1)); y=str2double(C(:,2)); %z=10;
% 
% CORR.x=A(:,15);  CORR.y=A(:,16);  CORR.z=A(:,17);
% SNR.x=A(:,11);  SNR.y=A(:,12);  SNR.z=A(:,13);
% corr=[CORR.x CORR.y CORR.z];
% snr=[SNR.x SNR.y SNR.z];
% A(corr<66)=NaN;
% A(snr<15)=NaN;
% A_out=fillmissing(A,"pchip");
% 
% %% Phase Space Threshold despiking (Goring & Nikora, 2002)
% 
% %% Turbulence
% %horizontal
% u=A_out(:,3);  v=A_out(:,4);  w=A_out(:,5);
% umean = mean(u);    vmean = mean(v);    wmean = mean(w);
% u_Tseries = u - umean;  v_Tseries = v - vmean; w_Tseries = w - wmean;
% 
% sample_variance_u = (length(u_Tseries)-1)*var(u_Tseries)/length(u_Tseries);
% sample_variance_v = (length(v_Tseries)-1)*var(v_Tseries)/length(v_Tseries); % estimation of population
% sample_variance_w = (length(w_Tseries)-1)*var(w_Tseries)/length(w_Tseries); 
% 
% urms=sqrt(sample_variance_u); vrms=sqrt(sample_variance_v); wrms=sqrt(sample_variance_w);
% 
% % cross moments (reynolds stresses)
% uwmean = mean(u_Tseries.*w_Tseries);
% RSS = abs(- 1000*uwmean);  % Reynolds shear stress
% uwcorr = uwmean/(sample_variance_w * sample_variance_u)^0.5;   % correlation coef
% 
% uvmean = mean(u_Tseries.*v_Tseries); vwmean = mean(v_Tseries.*w_Tseries);   
% TKE = 0.5*( sample_variance_u + sample_variance_v + sample_variance_w );
% %% higher order moments
% %TKE flux
% fku=0.5*(mean(u_Tseries.^3)+mean(u_Tseries .* v_Tseries.^2)+mean(u_Tseries .* w_Tseries.^2));
% fkv=0.5*(mean(v_Tseries .* u_Tseries.^2)+mean(v_Tseries.^3)+mean(v_Tseries .* w_Tseries.^2));
% fkw=0.5*(mean(u_Tseries.^2.*w_Tseries)+mean((w_Tseries .* v_Tseries.^2))+mean(w_Tseries.^3));
% 
% % 3rd order
% m30 = mean(u_Tseries.^3);
% m21 = mean(u_Tseries.^2 .* w_Tseries);
% m12 = mean(u_Tseries .* w_Tseries.^2);
% m03 = mean(w_Tseries.^3);
% 
% skewU = (length(u_Tseries)/(length(u_Tseries)-1))*mean(u_Tseries.^3)/sample_variance_u^(3/2);
% skewV = (length(v_Tseries)/(length(v_Tseries)-1))*mean(v_Tseries.^3)/sample_variance_v^(3/2);
% skewW = (length(w_Tseries)/(length(w_Tseries)-1))*mean(w_Tseries.^3)/sample_variance_w^(3/2);
% 
% % 4th order
% m40 = mean(u_Tseries.^4);
% m31 = mean(u_Tseries.^3 .* w_Tseries);
% m22 = mean(u_Tseries.^2 .* w_Tseries.^2);
% m13 = mean(u_Tseries .* w_Tseries.^3 );
% m04 = mean(w_Tseries.^4);
% 
% kurtU = (length(u_Tseries)/(length(u_Tseries)-1))*mean(u_Tseries.^4)/sample_variance_u^2 - 3;
% kurtV = (length(v_Tseries)/(length(v_Tseries)-1))*mean(v_Tseries.^4)/sample_variance_v^2 - 3;
% kurtW = (length(w_Tseries)/(length(w_Tseries)-1))*mean(w_Tseries.^4)/sample_variance_w^2 - 3;
% 
% %% Quadrants
% vel=[u_Tseries w_Tseries];
% H=0;
% threshold=H*abs(uwmean);
% ind=abs(u_Tseries.*w_Tseries)<threshold;
% vel(ind,:)=NaN;
% q1=find(vel(:,1)>0&vel(:,2)>0); q2=find(vel(:,1)<0&vel(:,2)>0);
% q3=find(vel(:,1)<0&vel(:,2)<0); q4=find(vel(:,1)>0&vel(:,2)<0);
% q=length(q1)+length(q2)+length(q3)+length(q4);
% RS1=length(q1)/q; RS2=length(q2)/q; RS3=length(q3)/q; RS4=length(q4)/q;
% %successive events time gap and muliply with 0.04=mean bursting period
% e_s=RS2/RS4;
% Te=(mean(diff(q2))*0.04);
% Ts=(mean(diff(q4))*0.04);
% % figure(1)
% % plot(u_Tseries, w_Tseries, 'r.')
% % hold on
% % plot(u_Tseries,abs(uwmean)./u_Tseries)
% %% Anisotropy
% RST=[mean(u_Tseries.^2)/(2*TKE) - 1/3, mean(u_Tseries.*v_Tseries)/(2*TKE), mean(u_Tseries.*w_Tseries)/(2*TKE) - 1/3
%     mean(u_Tseries.*v_Tseries)/(2*TKE), mean(v_Tseries.^2)/(2*TKE) - 1/3, mean(u_Tseries.*w_Tseries)/(2*TKE)
%     mean(u_Tseries.*w_Tseries)/(2*TKE), mean(w_Tseries.*v_Tseries)/(2*TKE), mean(w_Tseries.^2)/(2*TKE) - 1/3];
% eigen=eig(RST);
% I=eigen(1)+eigen(2)+eigen(3);    
% II=eigen(1)*eigen(2)+eigen(1)*eigen(3)+eigen(3)*eigen(2) ;     
% III=eigen(1)*eigen(2)*eigen(3);
% IF=1+27*III+9*II;
% Ad=(III/2)/(-II/3)^1.5;
tmp=eig(RST);
lambda = sort(tmp,'descend');     c = zeros(1,3);
c(1) = lambda(1)-lambda(2); c(2) = 2.*(lambda(2)-lambda(3)); c(3) = 3.*lambda(3)+1;
c=[0.2 0.4 0.4;0.1 0.6 0.3];
xc = [0.5, -0.5,           0];    yc = [   0,   0, sqrt(3)*0.5];
x1c=xc*c'; y1c=yc*c';z1c=c(:,1)';


x1= -0;
x2= 1;
x3=0.5;
y1=0;
y2=0;
y3=(sqrt(3)/2);
line([x1,x2], [y1, y2], 'Color', 'r');
hold on;
line([x2,x3], [y2, y3], 'Color', 'r');
line([x3,x1], [y3, y1], 'Color', 'r');
hold on;
axis off;
for i=1:length(y1c)
scatter(x1c(i),y1c(i),MarkerFaceColor=[c(i,1) c(i,2) c(i,3)])
end

for i=1:length(y1c)
contourf(x1c(i),y1c(i),z1c(i),cmap=[c(i,1) c(i,2) c(i,3)])
end
%eigen values
%eigen values
p=eig(b_mat);
s=b_mat^2;
inv_1=0;
inv_2=((p(1)*p(2))+(p(2)*p(3))+(p(3)*p(1)));
%inv_2=0.5*((trace(b_mat))^2-(trace(s)));
inv_3=det(b_mat);
eta=((inv_3/2)^(1/3));
nu=((-(inv_2)/3))^(0.5);
x1=linspace(0,-(1/6),1000);
y1=linspace(0,(1/6),1000);
plot(x1,y1);
hold on
x2=linspace(0,(1/3),1000);
y2=linspace(0,(1/3),1000);
plot(x2,y2);
hold on
x3=linspace(-(1/6),(1/3),1000);
y3=((1/27)+(2.*(x3).^3)).^0.5;
plot(x3,y3)
hold on 
grid on
plot(eta,nu,'r*');
title('LUMLEYS TRIANGLE')
i_f=1+27*inv_3+9*inv_2;
% load the data


% Set up a time variable
acq_freq = 25; %sampling frequency is 25Hz
dt = 1./acq_freq; %time interval between successive data points
t = 0:dt:(length(u)-1)*dt; %time of each sample

ubar = mean(u);
u = u - ubar;

%% 1. Calculate the energy spectra using an FFT. 

% Spectra (using the whole sample)
N = length(u); % number of samples
f = [0:N-1]./(N*dt); % determine the frequency range
uh = fft(u)/N; % do the FFT
eh = (conj(uh).*uh)*dt*N; % translate this into a spectrum (multiplying by N*dt ensures that the units are correct)
f = f(1:(N-1)/2); eh = eh(1:(N-1)/2); % cut off the data at the aliasing frequency
figure; loglog(f,eh); % plot the data on a log log scale
xlabel('frequency [Hz]')
ylabel('power spectral density')
% %% output
% k=[x y umean vmean wmean urms vrms wrms RSS uwcorr TKE fku fkv fkw m30 m21 m12 m03 skewU skewV skewW m40 m31 m22 m13 m04 kurtU kurtV kurtW RS1 RS2 RS3 RS4 Te Ts IF Ad];
% fileID = fopen("0.2_bed.dat","a+");
% fprintf(fileID,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n',k);
% end
% %% contours
% % N = 500;                           % Number Of Points Desired
% % xv = linspace(min(x), max(x), N);
% % yv = linspace(min(y), max(y), N);
% % [X,Y] = ndgrid(xv, yv);
% % Z = griddata(x, y, z, X, Y);
% % figure
% % contourf(X, Y, Z); shading interp;
% % hold on
% % rectangle('Position',[1 2 5 6])
% % hold on
% % rectangle('Position',[1 2 5 6])
% % axis('equal')
