% function fig_3_data
cc
% fig_3_data.m
%
% Generates data for figure 3 - note, due to the random number generator,
% the data will differ on each run. Uncomment line 75 to save data.

% Fixed parameters
N = 1e5;
I0 = 100;
R0 = 3;
gamma = 1/5;
alpha1 = gamma/49;
alpha2 = alpha1;
beta1 = R0*(alpha1+gamma)/N;
beta2 = beta1*1.5;
xi = 1e-5/gamma;
MaxTime = 365;
maxSteps = 1e7;
simtotal = 10;
NPIthreshold_on = 0;
NPIthreshold_off = 0;
c = 1;

% Variables
R = [0,0.5,0.75];

% Figure
figure(3)
clf
set(gcf,'color','w')

% Run simulations
for k=1:length(R)
    r = R(k);
    % Outputs
    I1_OUT = zeros(MaxTime+1,simtotal);
    I2_OUT = zeros(MaxTime+1,simtotal);
    I12_OUT = zeros(MaxTime+1,simtotal);
    I21_OUT = zeros(MaxTime+1,simtotal);
    C_OUT = zeros(MaxTime+1,simtotal);
    for s=1:simtotal
        [T, ~, I1, I2, C, I12, I21, ~, ~, ~, ~, ~] = variantevo(MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, c, r, N, I0, NPIthreshold_on, NPIthreshold_off, maxSteps, 1);
        
        for t=0:MaxTime
            t_temp = find(T>=t,1);
            if(~isempty(t_temp))
                I1_OUT(t+1,s) = I1(t_temp);
                I2_OUT(t+1,s) = I2(t_temp);
                I12_OUT(t+1,s) = I12(t_temp);
                I21_OUT(t+1,s) = I21(t_temp);
                C_OUT(t+1,s) = C(t_temp);
            end
        end
        
        subplot(2,length(R),k)
        hold on
        plot(T,(I1+I21+C)/N,'color',[0,0,1,0.25],'linewidth',1)
        plot(T,(I2+I12+C)/N,'color',[1,0,0,0.25],'linewidth',1)
        set(gca,'yscale','log')
        xlim([0,MaxTime])
        ylim([0.9/N,1])
        box on
        
        subplot(2,length(R),k+length(R))
        hold on
        plot(T,(I2+I12+C)./(I1+I2+I12+I21+C),'color',[0,0,0,0.25],'linewidth',1)
        set(gca,'yscale','log')
        xlim([0,MaxTime])
        ylim([0.9/N,1])
        box on
        set(gca,'ytick',logspace(log10(1/N),0,6))
    end
    
    % Uncomment line below to save data
    save(strcat('Data/fig_3_data_',num2str(r),'.mat'))
end

