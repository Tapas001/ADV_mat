clc
clear
%% Data preprocessing
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = [10, Inf];
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["TimeSeconds", "Position", "Flag", "Vx_0", "Vy_0", "Vz_0", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"];
opts.SelectedVariableNames = ["TimeSeconds", "Vx_0", "Vy_0", "Vz_0"];
opts.VariableTypes = ["double", "string", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];


% Import the data
file = "C:\experiments\T_head\x_15y_6\45_30_75.vf"; % change the file path according to you
[filepath,name,ext] = fileparts(file);
A = readmatrix(file,opts);
C = strsplit(name,'_');
x=str2double(C(:,1)); y=str2double(C(:,2)); %z=10;
A_out=fillmissing(A,"linear");
%% velocity Time series
u=A_out(:,2);  v=A_out(:,3);  w=A_out(:,4);
umean = mean(u);    vmean = mean(v);    wmean = mean(w);
u_Tseries = u - umean;  v_Tseries = v - vmean; w_Tseries = w - wmean;
%% Frequency calculations
%no of velocity samples
no_samples=size(u,1);
% Catergorizing velocities into different quadrants
q1=find(u_Tseries>=0 & w_Tseries >=0); q2=find(u_Tseries<0 & w_Tseries >=0);
q3=find(u_Tseries<=0 & w_Tseries<0); q4=find(u_Tseries>0 & w_Tseries<0);
N1=size(q1,1); N2=size(q2,1); N3=size(q3,1); N4=size(q4,1);
NQ=[N1; N2; N3; N4]; %number of points per quadrant
f = [N1/no_samples; N2/no_samples; N3/no_samples; N4/no_samples];
time_p = [N1*0.04; N2*0.04; N3*0.04; N4*0.04]; % For 25Hz = (1/25 = 0.04) ADV Data, May vary according to Frequency measurement
%% Stress Fraction Calculation
 urms=sqrt(mean(u_Tseries.*u_Tseries)); wrms=sqrt(mean(w_Tseries.*w_Tseries));
uw_denom =(sum(u_Tseries.* w_Tseries)) ; % Will be used as denominator for calculating Stress Fraction
threshold=urms*wrms; %threshold Value
H=0:1:8; % you can change Hole sizes
ukq1=zeros(size(u,1),size(H,2)); ukq2=zeros(size(u,1),size(H,2));
ukq3=zeros(size(u,1),size(H,2)); ukq4=zeros(size(u,1),size(H,2));
% 
 wkq1=zeros(size(u,1),size(H,2)); wkq2=zeros(size(u,1),size(H,2));
 wkq3=zeros(size(u,1),size(H,2)); wkq4=zeros(size(u,1),size(H,2));
% 
 for k=1:size(H,2)
     tauthresh=threshold;
     kappa=H(k);
     A1{k}=find(abs(u_Tseries.*w_Tseries)>(kappa*tauthresh) & w_Tseries>0 & u_Tseries>0);
     A2{k}=find(abs(u_Tseries.*w_Tseries)>(kappa*tauthresh) & w_Tseries>0 & u_Tseries<0);
     A3{k}=find(abs(u_Tseries.*w_Tseries)>(kappa*tauthresh) & w_Tseries<0 & u_Tseries<0);
     A4{k}=find(abs(u_Tseries.*w_Tseries)>(kappa*tauthresh) & w_Tseries<0 & u_Tseries>0);

% calculating stress fraction     
     S1(k)=sum((u_Tseries(A1{k}).*w_Tseries(A1{k})),1)/((uw_denom));
     S2(k)=sum((u_Tseries(A2{k}).*w_Tseries(A2{k})),1)/((uw_denom));
     S3(k)=sum((u_Tseries(A3{k}).*w_Tseries(A3{k})),1)/((uw_denom));
     S4(k)=sum((u_Tseries(A4{k}).*w_Tseries(A4{k})),1)/((uw_denom));
     n1(k)=size(A1{k},1); n2(k)=size(A2{k},1); n3(k)=size(A3{k},1); n4(k)=size(A4{k},1);
     
 end
 s1_0 = S1(1);
 s2_0 = S2(1);
 s3_0 = S3(1);
 s4_0 = S4(1);
 %% Storing values into dat file
% Remember the output will be stored only for H=0 
 k=[x y s1_0 s2_0 s3_0 s4_0 f(1,1) f(2,1) f(3,1) f(4,1) time_p(1,1) time_p(2,1) time_p(3,1) time_p(4,1)];
fileID = fopen("tapu1.dat","a+");
fprintf(fileID,'%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n',k);

 %% Plotting
    for i=1:length(A1)
        ufmaxi=max(abs(u_Tseries)); wfmaxi=max(abs(w_Tseries));
        xtik=ceil(ufmaxi*10)/10; ytik=ceil(wfmaxi*10)/10;
        figure, plot(u_Tseries, w_Tseries, 'w+'), hold on, plot(u_Tseries(A1{i}),w_Tseries(A1{i}), 'k+'), hold on,...
            plot(u_Tseries(A2{i}),w_Tseries(A2{i}), 'r+'),  hold on, plot(u_Tseries(A3{i}),w_Tseries(A3{i}), 'b+'), hold on,...
            plot(u_Tseries(A4{i}),w_Tseries(A4{i}), 'g+'),...
            axis([-ufmaxi ufmaxi -wfmaxi wfmaxi]);
               set(gca,'Fontsize',12, 'linewidth',2)
        title(['H=',num2str(kappa)],'Fontsize',12); %'
        xlabel('$u''$ $\mathrm{ms^{-1}}$','Fontsize',14,'Fontname','Times','Interpreter','Latex');
        ylabel('$v''$ $\mathrm{ms^{-1}}$','Fontsize',14,'Fontname','Times','Interpreter','Latex');
        set(gca,'Fontsize',11,'Fontname','Times'),
    end
    
ymax=max([max(max(abs(S1))) max(max(abs(S2))) max(max(abs(S3))) max(max(abs(S4)))]);
%ymax=1;%max(ymax);
xmax=max(H);
xmin=min(H);
ymin=0;
figure, subplot(222), plot(H, abs(S1), 'k+-','linewidth',2,'MarkerSize',10), axis([xmin xmax ymin ymax]),...
    set(gca,'Fontsize',12, 'linewidth',2), ...
    xlabel('$H$','Fontsize',14,'Fontname','Times','Interpreter','Latex'), ...
    ylabel('$S_{1}$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),...
    subplot(221), plot(H,abs(S2), 'ro-','linewidth',2,'MarkerSize',10), set(gca,'XDir','reverse'),axis([xmin xmax ymin ymax]),...
    set(gca,'Fontsize',12, 'linewidth',2),...
    xlabel('$H$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),...
    ylabel('$S_{2}$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),......
    subplot(223), plot(H,abs(S3), 'bx-','linewidth',2,'MarkerSize',10), set(gca,'YDir','reverse'), set(gca,'XDir','reverse'),...
    axis([xmin xmax ymin ymax]),set(gca,'Fontsize',12, 'linewidth',2),...
    xlabel('$H$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),...
    ylabel('$S_{3}$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),......
    subplot(224), plot(H,abs(S4), 'g*-','linewidth',2,'MarkerSize',10), set(gca,'YDir','reverse'),axis([xmin xmax ymin ymax]);
set(gca,'Fontsize',12, 'linewidth',2), set(gca,'Fontsize',12, 'linewidth',2),...
    xlabel('$H$','Fontsize',14,'Fontname','Times','Interpreter','Latex'),...
    ylabel('$S_{4}$','Fontsize',14,'Fontname','Times','Interpreter','Latex');

