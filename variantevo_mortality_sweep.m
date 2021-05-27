function variantevo_mortality_sweep

% variantevo_mortality_sweep.m
%
% Carries out the main parameter sweep for figures 4 and S1
% Dependency: variant_evo.c must be compiled first using "mex
% variant_evo.c"

% Fixed parameters
N = 1e5;
I0 = 100;
alpha1 = 0.005;
R0 = 3;
gamma = 1/5;
beta1 = R0*gamma/N;
xi = 1e-5/gamma;
MaxTime = 365;
maxSteps = 1e7;
fullOutput = 0;
simtotal = 1000;
C = 1;

% Variables
NPITHRESHOLD = [0,0.2];
R = linspace(0,1,101);
ALPHA2 = alpha1*[1,2,3,5];
BETA2MULT = 1.5;

% Carry out simulations
for m=1:length(NPITHRESHOLD)
    NPIthreshold = NPITHRESHOLD(m);
    for l=1:length(ALPHA2)
        alpha2 = ALPHA2(l);
        
        for s=1:length(BETA2MULT)
            beta2mult = BETA2MULT(s);
            beta2 = beta1*beta2mult;
            
            filename = strcat('Data/variantevo_mortality_',num2str(NPIthreshold),'_',num2str(beta2/beta1),'_',num2str(alpha2/alpha1),'.mat');
            
            if(exist(filename,'file'))
                disp('skipping:')
                disp(filename)
            else
                I2freqmax = zeros(length(C),length(R),simtotal);
                R1total = zeros(length(C),length(R),simtotal);
                R2total = zeros(length(C),length(R),simtotal);
                Rtotal = zeros(length(C),length(R),simtotal);
                Dtotal = zeros(length(C),length(R),simtotal);
                
                for j=1:length(R)
                    tic;
                    for i=1:length(C)
                        for k=1:simtotal
                            [~, ~, ~, ~, ~, ~, ~, R1total(i,j,k), R2total(i,j,k), Rtotal(i,j,k), Dtotal(i,j,k), I2freqmax(i,j,k)] = variantevo(MaxTime, alpha1, alpha2, beta1, beta2, gamma, xi, C(i), R(j), N, I0, NPIthreshold, maxSteps, fullOutput);
                        end
                    end
                    toc;
                    [i/length(C),j/length(R),l/length(ALPHA2),m/length(NPIthreshold),s/length(BETA2MULT)]
                end
                clear i j k
                
                save(filename)
            end
        end
    end
end