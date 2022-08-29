% Script for trimming the experimental backbone to be appropriate and
% saving in a convenient format for plotting. 

clear;

%% Load Data

% % Mode 1
% load('../EXPERIMENTAL_DATA/LTW_High_27Sept2021-2')
% Qmax = 13;
% Qmin = 0.3;
% export_name = 'LTW_Mode1Trimmed_27Sept2021';


% Mode 2
load('../EXPERIMENTAL_DATA/LTW_High_27Sept2021_Mode2')
Qmax = 42.3;
Qmin = 0.61;
export_name = 'LTW_Mode2Trimmed_27Sept2021';

%% Plot Data

figure; 

semilogx(AMP_avg, FRE_avg, 'k-'); hold on;
semilogx(AMP_avg, FRE_lower, 'k-');
semilogx(AMP_avg, FRE_upper, 'k-');
hold on;
xlabel('Amplitude');
ylabel('Frequency');

for ii = 1:length(Amps)
    
    semilogx(Amps{ii}, Freqs{ii}, '--');
end


figure; 

loglog(AMP_avg, DAM_avg, 'k-'); hold on;
loglog(AMP_avg, DAM_lower, 'k-');
loglog(AMP_avg, DAM_upper, 'k-');
hold on;
xlabel('Amplitude');
ylabel('Damping');

for ii = 1:length(Amps)
    
    semilogx(Amps{ii}, Damps{ii}, '--');
end


%% Trim Data


Qaccel = AMP_avg(AMP_avg > Qmin & AMP_avg < Qmax);

W_all = zeros(length(Qaccel), length(Amps));
Z_all = zeros(length(Qaccel), length(Amps));

for ii = 1:length(Amps)
    
    W_all(:, ii) = interp1(Amps{ii}, Freqs{ii}, Qaccel);
    Z_all(:, ii) = interp1(Amps{ii}, Damps{ii}, Qaccel);
    
end

BB.Qaccel = Qaccel;
BB.Wmean = mean(W_all, 2, 'omitnan');
BB.Zmean = mean(Z_all, 2, 'omitnan');

BB.Wmax = max(W_all, [], 2, 'omitnan');
BB.Zmax = max(Z_all, [], 2, 'omitnan');

BB.Wmin = min(W_all, [], 2, 'omitnan');
BB.Zmin = min(Z_all, [], 2, 'omitnan');

figure;
subplot(2,1,1)
semilogx(BB.Qaccel, BB.Wmean, 'k-'); hold on;
semilogx(BB.Qaccel, BB.Wmin, 'k-');
semilogx(BB.Qaccel, BB.Wmax, 'k-');

subplot(2,1,2)
loglog(BB.Qaccel, BB.Zmean, 'k-'); hold on;
loglog(BB.Qaccel, BB.Zmin, 'k-');
loglog(BB.Qaccel, BB.Zmax, 'k-');

% save(export_name, 'BB');



