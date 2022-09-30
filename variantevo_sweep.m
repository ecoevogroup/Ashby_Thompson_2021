function variantevo_sweep

% variantevo_sweep.m
%
% Carries out the main parameter sweep for figures 2 and 5 (and figures S1,
% S2, S4 and S5 by changing ALPHA1).
%
% Dependency: variant_evo.c must be compiled first using "mex
% variant_evo.c"

% Fixed parameters
N = 1e5;
I0 = 100;
R0 = 3;
gamma = 1/5;
xi = 1e-5/gamma;
MaxTime = 365;
maxSteps = 1e7;
fullOutput = 0;
simtotal = 1000;

% Variables
ALPHA1 = gamma./49;
NPITHRESHOLD_ON = [0,0.01];
NPITHRESHOLD_OFF = [0,0.002];
C = linspace(0,1,101);
R = linspace(0,1,101);
BETA2MULT = [0.5,1,1.5,2];

% Carry out simulations
for a=1:length(ALPHA1)
    alpha1 = ALPHA1(a);
    ALPHA2 = alpha1;
    beta1 = R0*(alpha1+gamma)/N;
    for m=1:length(NPITHRESHOLD_ON)
        NPIthreshold_on = NPITHRESHOLD_ON(m);
        NPIthreshold_off = NPITHRESHOLD_OFF(m);
        for l=1:length(ALPHA2)
            alpha2 = ALPHA2(l);
            
            for s=1:length(BETA2MULT)
                beta2mult = BETA2MULT(s);
                beta2 = beta1*beta2mult;
                
                filename = strcat('Data/variantevo_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat');
                
                if(exist(filename,'file'))
                    disp('skipping:')
                    disp(filename)
                else
                    I2freqmax = zeros(length(C),length(R),simtotal);
                    R1total = zeros(length(C),length(R),simtotal);
                    R2total = zeros(length(C),length(R),simtotal);
                    Rtotal = zeros(length(C),length(R),simtotal);
                    Dtotal = zeros(length(C),length(R),simtotal);
                    
                    save(filename)
                    
                    for j=1:length(R)
                        tic;
                        for i=1:length(C)
                            parfor k=1:simtotal
                                [~, ~, ~, ~, ~, ~, ~, R1total(i,j,k), R2total(i,j,k), Rtotal(i,j,k), Dtotal(i,j,k), I2freqmax(i,j,k)] = variantevo(MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, C(i), R(j), N, I0, NPIthreshold_on, NPIthreshold_off, maxSteps, fullOutput);
                            end
                        end
                        toc;
                        [i/length(C),j/length(R),a/length(ALPHA1), l/length(ALPHA2),m/length(NPITHRESHOLD_ON),s/length(BETA2MULT)]
                    end
                    clear i j k
                    
                    save(filename)
                end
            end
        end
    end
end

