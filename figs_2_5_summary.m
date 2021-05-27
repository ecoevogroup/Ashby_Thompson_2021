function figs_2_5_summary

% figs_2_5_summary.m
%
% Generates summary data for figures 2 and 5

% Fixed parameters
N = 1e5;
alpha1 = 0.005;
alpha2 = alpha1;
R0 = 3;
gamma = 1/5;
beta1 = R0*gamma/N;

% Variables
NPITHRESHOLD0 = [0,0.2];
C = linspace(0,1,101);
BETA2MULT00 = [0.75,1,1.5];

% Loop over each figure
for m=1:length(NPITHRESHOLD0)
    NPIthreshold = NPITHRESHOLD0(m);    
    
    % Loop over variation in the relative transmission rate
    for s0=1:length(BETA2MULT00)
        beta2mult = BETA2MULT00(s0);
        beta2 = beta1*beta2mult;
        
        filename = strcat('Data/variantevo_',num2str(NPIthreshold),'_',num2str(beta2/beta1),'_',num2str(alpha2/alpha1),'.mat');
        load(filename)

        median_variant_recovered = median(R2total+Rtotal,3)/N;
        median_deaths = median(Dtotal,3);        
        variant_emergence = mean(I2freqmax>0.1,3);
        
        clear R1total R2total Rtotal Dtotal I2freqmax filename beta2mult m s0
        save(strcat('Data/variantevo_summary_',num2str(NPIthreshold),'_',num2str(beta2/beta1),'_',num2str(alpha2/alpha1),'.mat'))
    end
end