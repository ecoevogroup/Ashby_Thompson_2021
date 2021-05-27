function figs_4_S1

% figs_4_S1.m
%
% Generates figures 4 and S1

% Fixed parameters
N = 1e5;
alpha1 = 0.005;
R0 = 3;
gamma = 1/5;
beta1 = R0*gamma/N;
beta2mult = 1.5;
beta2 = beta1*beta2mult;

% Variables
R = linspace(0,1,101);
ALPHA20 = alpha1*(1:3);
NPITHRESHOLD = [0,0.2];

% Figure labels
labs1 = {'(a)','(b)','(c)'};

% Loop over each figure
for k=1:2
    NPIthreshold = NPITHRESHOLD(k);
    if(k==1)
        figure(4)
    else
        figure(6)
    end
    
    % Setup figure
    clf
    set(gcf,'color','w')
    set(gcf,'PaperUnits','centimeters')
    xSize = 14; ySize = 4;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 100 xSize*50 ySize*50])
    
    % Loop over variation in the relative mortality rate
    for l=1:length(ALPHA20)
        alpha2 = ALPHA20(l);
        
        subplot(1,length(ALPHA20),l)
        hold on
        
        filename = strcat('Data/variantevo_mortality_',num2str(NPIthreshold),'_',num2str(beta2/beta1),'_',num2str(alpha2/alpha1),'.mat');
        load(filename)
        
        % Plot data
        OUT = prctile(Dtotal(1,:,:),[25,75],3);
        patch([R fliplr(R)],[max(0,OUT(:,:,1)) fliplr(OUT(:,:,2))],'k','facealpha',0.3,'linestyle','none');
        plot(R,median(Dtotal(end,:,:),3),'color','k','linewidth',2)
        
        box on
        ylim([0,900])
        text(0,max(get(gca,'ylim'))*1.1,labs1{l},'interpreter','latex','fontsize',14)
        if(l==2)
            xlabel('Strength of NPIs, $r$','interpreter','latex','fontsize',16)
        end
        if(l==1)
            ylabel('Total deaths (per 100k)','interpreter','latex','fontsize',16)
        end
        title(strcat('$\frac{\alpha_v}{\alpha_w}=',num2str(alpha2/alpha1),'$'),'interpreter','latex','fontsize',16)
    end
    
%     % Save figures
%     if(k==1)
%         save2pdf('fig_4.pdf')
%     else
%         save2pdf('fig_S1.pdf')
%     end
end