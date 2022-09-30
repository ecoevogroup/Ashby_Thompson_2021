function figs_4_6

% figs_4_6.m
%
% Generates figures 4 and 6

% Fixed parameters
N = 1e5;
R0 = 3;
gamma = 1/5;
alpha1 = gamma./49;
beta2mult = 2;
beta1 = R0*(alpha1+gamma)/N;
beta2 = beta1*beta2mult;

% Variables
ALPHA20 = alpha1*[1,2,3];
NPITHRESHOLD_ON = [0,0.01];
NPITHRESHOLD_OFF = [0,0.002];
R = linspace(0,1,101);

% Figure labels
labs1 = {'(a)','(b)','(c)'};

% Loop over each figure
for k=1:2
    NPIthreshold_on = NPITHRESHOLD_ON(k);
    NPIthreshold_off = NPITHRESHOLD_OFF(k);
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
    
    ymax = 0;
    % Loop over variation in the relative mortality rate
    for l=1:length(ALPHA20)
        alpha2 = ALPHA20(l);
        
        subplot(1,length(ALPHA20),l)
        hold on
        
        filename = strcat('Data/variantevo_mortality_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat');
        load(filename)
        
        % Plot data
        OUT = prctile(Dtotal(1,:,:),[25,75],3);
        patch([R fliplr(R)],[max(0,OUT(:,:,1)) fliplr(OUT(:,:,2))],'k','facealpha',0.3,'linestyle','none');
        plot(R,median(Dtotal(end,:,:),3),'color','k','linewidth',2)
        
        ymax = max(ymax,max(OUT(:,:,2)));
        
        box on
        if(l==2)
            xlabel('Strength of NPIs, $r$','interpreter','latex','fontsize',16)
        end
        if(l==1)
            ylabel('Total deaths (per 100k)','interpreter','latex','fontsize',16)
        end
        if(alpha1==alpha2)
            title(strcat('$\alpha_v=\alpha_w$'),'interpreter','latex','fontsize',16)
        else
            title(strcat('$\alpha_v=',num2str(alpha2/alpha1),'\alpha_w$'),'interpreter','latex','fontsize',16)
        end
    end
    
    ymax = 100*ceil(ymax/100);
    for l=1:length(ALPHA20)
        subplot(1,length(ALPHA20),l)
        ylim([0,ymax]);
        text(0,max(get(gca,'ylim'))*1.1,labs1{l},'interpreter','latex','fontsize',14)
    end
    
    %     % Save figures
    %     if(k==1)
    %         save2pdf('fig_4.pdf')
    %     else
    %         save2pdf('fig_6.pdf')
    %     end
end