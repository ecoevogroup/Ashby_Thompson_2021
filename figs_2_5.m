function figs_2_5

% figs_2_5.m
%
% Generates figures 2 and 5

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

% Figure labels
labs1 = {'(a)','(b)','(c)'};
labs2 = {'(d)','(e)','(f)'};
labs3 = {'(g)','(h)','(i)'};

% Loop over each figure
for m=1:length(NPITHRESHOLD0)
    NPIthreshold = NPITHRESHOLD0(m);    
    if(m==1)
        figure(2)
    else
        figure(5)
    end
    
    % Setup figure
    clf
    set(gcf,'color','w')
    set(gcf,'PaperUnits','centimeters')
    xSize = 10; ySize = 10;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 100 xSize*50 ySize*50])
    
    % Loop over variation in the relative transmission rate
    for s0=1:length(BETA2MULT00)
        beta2mult = BETA2MULT00(s0);
        beta2 = beta1*beta2mult;
        
        filename = strcat('Data/variantevo_summary_',num2str(NPIthreshold),'_',num2str(beta2/beta1),'_',num2str(alpha2/alpha1),'.mat');
        load(filename)
        
        % Plot data
        subplot(3,length(BETA2MULT00),s0)
        imagesc(median_variant_recovered)
%         imagesc(median(R2total+Rtotal,3)/N)
        set(gca,'ydir','normal')
        set(gca,'xtick',linspace(1,length(C),6),'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'ytick',linspace(1,length(C),6),'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'clim',[0,1])
        set(gca,'fontsize',10)
        temp = get(gca,'position');
        temp(1) = temp(1) - 0.05;
        set(gca,'position',temp);
        if(s0==length(BETA2MULT00))
            c=colorbar;
            temp = get(c,'position');
            temp(1) = temp(1) + 0.07;
            set(c,'position',temp);
            ylabel(c,{'Proportion infected','by variant'},'interpreter','latex','fontsize',14)
        end
        text(0,length(C)*1.1,labs1{s0},'interpreter','latex','fontsize',14)
        title(strcat('$\frac{\beta_v}{\beta_w}=',num2str(beta2mult),'$'),'interpreter','latex','fontsize',16)
        
        subplot(3,length(BETA2MULT00),s0+length(BETA2MULT00))
        imagesc(median_deaths)
%         imagesc(median(Dtotal,3))
        set(gca,'ydir','normal')
        set(gca,'xtick',linspace(1,length(C),6),'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'ytick',linspace(1,length(C),6),'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'clim',[0,900])
        set(gca,'fontsize',10)
        temp = get(gca,'position');
        temp(1) = temp(1) - 0.05;
        set(gca,'position',temp);
        if(s0==1)
            ylabel('Strength of cross immunity, $c$','interpreter','latex','fontsize',16)
        end
        if(s0==length(BETA2MULT00))
            c=colorbar;
            temp = get(c,'position');
            temp(1) = temp(1) + 0.07;
            set(c,'position',temp);
            ylabel(c,{'Total deaths','(per 100k)'},'interpreter','latex','fontsize',14')
        end
        text(0,length(C)*1.1,labs2{s0},'interpreter','latex','fontsize',14)
        
        subplot(3,length(BETA2MULT00),s0+2*length(BETA2MULT00))
        imagesc(variant_emergence)
        set(gca,'ydir','normal')
        set(gca,'xtick',linspace(1,length(C),6),'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'ytick',linspace(1,length(C),6),'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
        set(gca,'clim',[0,1])
        set(gca,'fontsize',10)
        temp = get(gca,'position');
        temp(1) = temp(1) - 0.05;
        set(gca,'position',temp);
        if(s0==2)
            xlabel('Strength of NPIs, $r$','interpreter','latex','fontsize',16)
        end
        if(s0==length(BETA2MULT00))
            c=colorbar;
            temp = get(c,'position');
            temp(1) = temp(1) + 0.07;
            set(c,'position',temp);
            ylabel(c,{'Variant emergence','probability'},'interpreter','latex','fontsize',14)
        end
        text(0,length(C)*1.1,labs3{s0},'interpreter','latex','fontsize',14) 
    end
    
%     % Save figures
%     if(m==1)
%         save2pdf('fig_2.pdf');
%     else
%         save2pdf('fig_5.pdf');
%     end
end