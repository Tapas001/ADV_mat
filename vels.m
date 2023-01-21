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
file = "c:\experiments\T_head\x_15y_6\45_30_77.Vf";
[filepath,name,ext] = fileparts(file);
A = readmatrix(file,opts);
C = strsplit(name,'_');
x=str2double(C(:,1)); y=str2double(C(:,2)); %z=10;
A_out=fillmissing(A,"linear");
%% getting data
%horizontal
u=A_out(:,2);  v=A_out(:,3);  w=A_out(:,4);
%% spectral analysis
[t f] = pspectrum(u,25);
%% different slopes
figure;
subplot(2,2,1)
loglog(f,t)
xlabel('Freqency');
ylabel('S(uu)');
title('Subplot1 :Energyspectrum');
subplot(2,2,2)
loglog(f,t)
y1 = [f.^(-(5/3))];
hold on
plot(f,y1)
xlabel('Freqency');
ylabel('S(uu)');
title('Subplot2: -5/3 slope');
subplot(2,2,3)
loglog(f,t)
y2 = 0.2*[f.^(-(3))];
hold on
plot(f,y2)
xlabel('Freqency');
ylabel('S(uu)');
title('Subplot3:-3 slope')
subplot(2,2,4);
loglog(f,t)
y3 = 0.6*[f.^(-(1))];
hold on
plot(f,y3)
xlabel('Freqency');
ylabel('S(uu)');
title('Subplot4 :-1 slope');
% premultiplied spectra
s = t.*f;
% figure;
% semilogx(f,s)
% xlabel('Freqency');
% ylabel('S(uu)* Frequency');
% title('Premultipliedspectra');
%% ADV data validation
r = sum(s);
sum = 0;
L = zeros(length(s),1);
for i = 1:length(s)
    sum = sum + s(i);
    L(i) = sum;
end
p = L./r;
figure;
subplot(2,2,1)
plot(f,p)
xlabel('Freqency');
ylabel('Gamma');
title('Subplot1: Energyfraction normal scale');
subplot(2,2,2)
semilogx(f,p)
xlabel('Freqency');
ylabel('Gamma');
title('Subplot2:Energy fraction semilog scale')
subplot(2,2,[3,4]);
semilogx(f,s)
xlabel('Freqency');
ylabel('S(uu)* Frequency');
title('Subplot3 : Premultipliedspectra');

