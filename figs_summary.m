function figs_summary

% figs_summary.m
%
% Generates summary data for figures 2, 5, S1-2, S4-5

% Fixed parameters
N = 1e5;
R0 = 3;
gamma = 1/5;

% Variables
ALPHA10 = gamma./[499,49,4];
NPITHRESHOLD_ON0 = [0,0.01];
NPITHRESHOLD_OFF0 = [0,0.002];
C = linspace(0,1,101);
BETA2MULT00 = [0.5,1,1.5,2];

% Loop over each figure
for a=1:length(ALPHA10)
    alpha1 = ALPHA10(a);
    alpha2 = alpha1;
    beta1 = R0*(alpha1+gamma)/N;
    for m=1:length(NPITHRESHOLD_ON0)
        NPIthreshold_on = NPITHRESHOLD_ON0(m);
        NPIthreshold_off = NPITHRESHOLD_OFF0(m);
        
        % Loop over variation in the relative transmission rate
        for s0=1:length(BETA2MULT00)
            beta2mult = BETA2MULT00(s0);
            beta2 = beta1*beta2mult;
            
            filename = strcat('Data/variantevo_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat');
            if(exist(filename,'file'))
                load(filename)
                
                median_variant_recovered = median(R2total+Rtotal,3)/N;
                median_deaths = median(Dtotal,3);
                variant_emergence = mean(I2freqmax>0.1,3);
                
                clear R1total R2total Rtotal Dtotal I2freqmax filename beta2mult m s0 a BETA2MULT NPITHRESHOLD_ON NPITHRESHOLD_OFF
                save(strcat('Data/variantevo_summary_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat'))
            end
        end
    end
end