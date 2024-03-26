clc;
clear all


%% 
for z = 1:2:104
     A = z; % select fluorescent data column for LCGR
     B = z+1; % select fluorescent data column for TXR
  
% A=11;
% B=12;

%% enter melt prePCR data 

% enter file name for data file
%%%data = readmatrix('time format _ formatted raw tab only 228 pt dont use button just paste into tab 1 --- MACRO 2023-10-04 Expt 236 pcr and post melt cont 025 ramp');
data = readmatrix('E243 temp formatted raw tab only');

% this is a macro formatted quant studio 5 melt curve result data file!!!!!

% note melts have 228 points using 0.025C/s continuous acq 
%%%% this melt has 226 points in E243!!!

% ******************************************************************
% $$$$$$ ENTER NUMBER OF DATA COLUMNS IN FIRST ENTRY FOR `42` $$$$$ 
% ******************************************************************

data = data(:,1:208);
x = data(1:length(data),1:2:end); % organize temp(C) data
y = data(1:length(data),2:2:end); % organize raw fluorescence data

%%


%% Li Griffith Univ Approach
k = 2; % The polynomial order k is set as 2

RecWindow = 5; % RecWindow is the number of samples within 1 C. 
%https://www.researchgate.net/profile/Cameron-Gundry/publication/10743956_High-Resolution_Genotyping_by_Amplicon_Melting_Analysis_Using_LCGreen/links/09e415085b75752b21000000/High-Resolution-Genotyping-by-Amplicon-Melting-Analysis-Using-LCGreen.pdf
%wittwer: Derivative melting curve plots were calculated from
% the Savitsky–Golay polynomials at each point (17).
% Savitsky–Golay analysis used a second-degree polynomial and a data window including all points within a 1 °C
% interval.

% Calculated from 228pts/1200s = 0.19pts/s 
% 5.24s/0.131C = 40s/C --- 40s/C*0.19pts/s = 7.6pts/C --- round down to 7
% for integer frame length

f = (RecWindow)*2 + 1; % f is frame size
[b,g] = sgolay(k,f); % Savitzky–Golay FIR smoothing filter b and also the matrix g of differentiation filters

F_LCGR = y(:,A); 
F_TXR = y(:,B);

% an array of the measured raw fluorescence melting curve data, F()
% The n-th data point at the negative first derivative curve is found by 
% apply the MATLAB dot() function to return the scalar product of the two vector segments
%DerivFluo(n) = -dot(g(:,2),F(n-RecWindow : n+RecWindow));

for n=1+RecWindow:length(F_LCGR)-RecWindow
    DerivFluo_LCGR_rows(n-RecWindow) = -dot(g(:,2),F_LCGR(n-RecWindow : n+RecWindow));
end 
DerivFluo_LCGR_columns = DerivFluo_LCGR_rows';

for n=1+RecWindow:length(F_TXR)-RecWindow
    DerivFluo_TXR_rows(n-RecWindow) = -dot(g(:,2),F_TXR(n-RecWindow : n+RecWindow));
end 
DerivFluo_TXR_columns = DerivFluo_TXR_rows';
%%

% find Tms 
[~, Tm_LCGR_timelocation] = max(DerivFluo_LCGR_columns(1:212, :));
Tm_LCGR_timelocation = Tm_LCGR_timelocation + RecWindow;
[~, Tm_TXR_timelocation] = min(DerivFluo_TXR_columns(1:212, :));
Tm_TXR_timelocation = Tm_TXR_timelocation + RecWindow;

% print Tms
Tm_LCGR = x(Tm_LCGR_timelocation+0);
Tm_TXR = x(Tm_TXR_timelocation+0);

DerivFluo_columns_summary(:,A) = DerivFluo_LCGR_columns;
DerivFluo_columns_summary(:,B) = DerivFluo_TXR_columns;


Tm_summary(:,A) = Tm_LCGR;
Tm_summary(:,B) = Tm_TXR;
x_summary(:,A) = x(1+RecWindow:length(x)-RecWindow,A);
x_summary(:,B) = x(1+RecWindow:length(x)-RecWindow,B);
% Tm_dif_summary (:,:) = Tm_LCGR-Tm_TXR;

%%
% figure
%%% plot raw curve
    % hold on
    % plot(x(:,A),F_LCGR, 'g', x(:,B),F_TXR*-18, 'r')
    % yyaxis right

%% plot deriv curve
% plot(x(1+RecWindow:length(x)-RecWindow,A),DerivFluo_LCGR_columns, 'g', x(1+RecWindow:length(x)-RecWindow,B),DerivFluo_TXR_columns*-18, 'r')
% % xline(Tm_LCGR(A),'g')
% % xline(Tm_TXR(B),'r')
% xline(Tm_LCGR,'g')
% xline(Tm_TXR,'r')

%% 
%  plot(x,y); % plot raw data

%%
%  x2 = x(1:length(x)-1,:); % make x space with one less number for use in derivative

% %% plot LCGR+ including Tm's for smooth vs unedited analysis approaches - broken down samples 4,5,6
% 
% figure; % create figure with subplots
% 
% % subplot of sample 4 WT 10^6
% subplot(3,1,1); % 2 rows, 1 column, first subplot
% plot(x2(:,4),dydx_unedited_lowpass(:,4),'g',x2(:,4),dydx_lowpass(:,4),'r'); %plot only sample 4 derivs 
% xline(Tm_premelt_unedited(4),'g',Tm_premelt_unedited(4))
% xline(Tm_premelt_smooth(4),'r',Tm_premelt_smooth(4))
% xline(82.392,'k') % quant output Tm for sample 4 in G1
% title('Sample 4 WT Rep 1 @ 10^6: Red (Sgolay Smoothing on Raw + Low Pass Filter Derivative) vs. Green (Low Pass Filter Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of sample 5 WT 10^6
% subplot(3,1,2); % 2 rows, 1 column, first subplot
% plot(x2(:,5),dydx_unedited_lowpass(:,5),'g',x2(:,5),dydx_lowpass(:,5),'r'); %plot only sample 5 derivs 
% xline(Tm_premelt_unedited(5),'g',Tm_premelt_unedited(5))
% xline(Tm_premelt_smooth(5),'r',Tm_premelt_smooth(5))
% xline(82.392,'k') % quant output Tm for sample 5 in G2
% title('Sample 5 WT Rep 2 @ 10^6: Red (Sgolay Smoothing on Raw + Low Pass Filter Derivative) vs. Green (Low Pass Filter Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of sample 6 WT 10^6
% subplot(3,1,3); % 2 rows, 1 column, first subplot
% plot(x2(:,6),dydx_unedited_lowpass(:,6),'g',x2(:,6),dydx_lowpass(:,6),'r'); %plot only sample 6 derivs 
% xline(Tm_premelt_unedited(6),'g',Tm_premelt_unedited(6))
% xline(Tm_premelt_smooth(6),'r',Tm_premelt_smooth(6))
% xline(82.509,'k') % quant output Tm for sample 6 in G3
% % title('Sample 6 WT Rep 2 @ 10^6: Red (Sgolay Smoothing on Raw + Low Pass Filter Derivative) vs. Green (Low Pass Filter Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');

% %% plot TXR vs LCGR+ without Tm's for unedited lowpass analysis on rep 1's per target type
% 
% figure; % create figure with subplots
% 
% % LCGR+ in green and TXR in red
% 
% % subplot of F1 MEP183
% subplot(5,2,1); % 5 rows, 2 columns, first subplot
% plot(x2(:,7),dydx_unedited_lowpass(:,7),'g',x2(:,8),dydx_unedited_lowpass(:,8)*18,'r'); 
% xline(Tm_premelt_unedited(7),'g')
% xline(Tm_premelt_unedited(8),'r')
% xline(82.5031433105469,'k') % QS5 output Tm for F1 LCGR+
% title('WT (MEP183, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of F4 MEP184
% subplot(5,2,2); 
% plot(x2(:,13),dydx_unedited_lowpass(:,13),'g',x2(:,14),dydx_unedited_lowpass(:,14)*18,'r'); 
% xline(Tm_premelt_unedited(13),'g')
% xline(Tm_premelt_unedited(14),'r')
% xline(82.4570465087891,'k') % QS5 output Tm for F4 LCGR+
% title('S315T (MEP184, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of F7 MEP185
% subplot(5,2,3); 
% plot(x2(:,19),dydx_unedited_lowpass(:,19),'g',x2(:,20),dydx_unedited_lowpass(:,20)*18,'r'); 
% xline(Tm_premelt_unedited(19),'g')
% xline(Tm_premelt_unedited(20),'r')
% xline(81.8079833984375,'k') % QS5 output Tm for F7 LCGR+
% title('S315N (MEP185, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of F10 MEP186
% subplot(5,2,4); 
% plot(x2(:,25),dydx_unedited_lowpass(:,25),'g',x2(:,26),dydx_unedited_lowpass(:,26)*18,'r'); 
% xline(Tm_premelt_unedited(25),'g')
% xline(Tm_premelt_unedited(26),'r')
% xline(81.7529067993164,'k') % QS5 output Tm for F10 LCGR+
% title('S315I (MEP186, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of G1 MEP187
% subplot(5,2,5); 
% plot(x2(:,31),dydx_unedited_lowpass(:,31),'g',x2(:,32),dydx_unedited_lowpass(:,32)*18,'r'); 
% xline(Tm_premelt_unedited(31),'g')
% xline(Tm_premelt_unedited(32),'r')
% xline(81.5819396972656,'k') % QS5 output Tm for G1 LCGR+
% title('S315R (MEP187, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of G4 MEP188
% subplot(5,2,6); 
% plot(x2(:,37),dydx_unedited_lowpass(:,37),'g',x2(:,38),dydx_unedited_lowpass(:,38)*18,'r'); 
% xline(Tm_premelt_unedited(37),'g')
% xline(Tm_premelt_unedited(38),'r')
% xline(83.2475128173828,'k') % QS5 output Tm for G4 LCGR+
% title('S315G (MEP188, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of G7 MEP189
% subplot(5,2,7); 
% plot(x2(:,43),dydx_unedited_lowpass(:,43),'g',x2(:,44),dydx_unedited_lowpass(:,44)*18,'r'); 
% xline(Tm_premelt_unedited(43),'g')
% xline(Tm_premelt_unedited(44),'r')
% xline(82.1362380981445,'k') % QS5 output Tm for G7 LCGR+
% title('S315L (MEP189, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of G10 MEP197
% subplot(5,2,8); 
% plot(x2(:,49),dydx_unedited_lowpass(:,49),'g',x2(:,50),dydx_unedited_lowpass(:,50)*18,'r'); 
% xline(Tm_premelt_unedited(49),'g')
% xline(Tm_premelt_unedited(50),'r')
% xline(81.3587493896484,'k') % QS5 output Tm for G10 LCGR+
% title('S315T + A312V (MEP197, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of H1 MEP198
% subplot(5,2,9); 
% plot(x2(:,55),dydx_unedited_lowpass(:,55),'g',x2(:,56),dydx_unedited_lowpass(:,56)*18,'r'); 
% xline(Tm_premelt_unedited(55),'g')
% xline(Tm_premelt_unedited(56),'r')
% xline(81.6477432250977,'k') % QS5 output Tm for H1 LCGR+
% title(' S315T + G316D (MEP198, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');
% 
% % subplot of H4 MEP199
% subplot(5,2,10); 
% plot(x2(:,61),dydx_unedited_lowpass(:,61),'g',x2(:,62),dydx_unedited_lowpass(:,62)*18,'r'); 
% xline(Tm_premelt_unedited(61),'g')
% xline(Tm_premelt_unedited(62),'r')
% xline(80.9419937133789,'k') % QS5 output Tm for H4 LCGR+
% title(' S315T + A312V + G316D (MEP199, Rep 1): LCGR+ vs TXR(scaled*18) via Low Pass Filter on Derivative');
% xlabel('Temperature (ºC)');
% ylabel('abs(dy/dx)');

%% 

% % plot TXR vs LCGR+ without Tm's for unedited lowpass analysis on rep 1's per target type
% 
% figure; % create figure with subplots
% 
% % LCGR+ in green and TXR in red
% 
% % subplot of F1 MEP183
% subplot(3,1,1); % 5 rows, 2 columns, first subplot
% plot(x2(:,7),dydx_unedited_lowpass(:,7),'g',x2(:,8),dydx_unedited_lowpass(:,8)*18,'r'); 
% xline(Tm_premelt_unedited(7),'g')
% xline(Tm_premelt_unedited(8),'r')
% %xline(82.5031433105469,'k') % QS5 output Tm for F1 LCGR+
% title('WT (MEP183, Rep 1): LCGR+ vs TXR(scaled*18)');
% xlabel('Time (sec)');
% ylabel('abs(dy/dx)');
% % xlim([600 800]);
% % ylim([0 7000]);
% 
% % subplot of F4 MEP184
% subplot(3,1,2); 
% plot(x2(:,13),dydx_unedited_lowpass(:,13),'g',x2(:,14),dydx_unedited_lowpass(:,14)*18,'r'); 
% xline(Tm_premelt_unedited(13),'g')
% xline(Tm_premelt_unedited(14),'r')
% %xline(82.4570465087891,'k') % QS5 output Tm for F4 LCGR+
% title('S315T (MEP184, Rep 1): LCGR+ vs TXR(scaled*18)');
% xlabel('Time (sec)');
% ylabel('abs(dy/dx)');
% % xlim([600 800]);
% % ylim([0 7000]);
% 
% % subplot of H4 MEP199
% subplot(3,1,3); 
% plot(x2(:,61),dydx_unedited_lowpass(:,61),'g',x2(:,62),dydx_unedited_lowpass(:,62)*18,'r'); 
% xline(Tm_premelt_unedited(61),'g')
% xline(Tm_premelt_unedited(62),'r')
% %xline(80.9419937133789,'k') % QS5 output Tm for H4 LCGR+
% title(' S315T + A312V + G316D (MEP199, Rep 1): LCGR+ vs TXR(scaled*18)');
% xlabel('Time (sec)');
% ylabel('abs(dy/dx)');
% % xlim([600 800]);
% % ylim([0 7000]);
% 

%%
 end

%  %%
% figure
% plot(x(1+RecWindow:length(x)-RecWindow,7),DerivFluo_LCGR_columns, 'g', x(1+RecWindow:length(x)-RecWindow,8),DerivFluo_TXR_columns*-18, 'r')
% xline(Tm_LCGR(7),'g')
% xline(Tm_TXR(8),'r')