function figs_S1_S2_S4_S5

% figs_S1_S2_S4_S5.m
%
% Generates figures S1, S2, S4 and S5

% Fixed parameters
N = 1e5;
R0 = 3;
gamma = 1/5;

% Variables
ALPHA100 = [gamma/499,gamma/4];
NPITHRESHOLD_ON0 = [0,0.01];
NPITHRESHOLD_OFF0 = [0,0.002];
C = linspace(0,1,101);
BETA2MULT = [0.5,1,1.5,2];

% Figure labels
labs1 = {'(a)','(b)','(c)','(d)'};
labs2 = {'(e)','(f)','(g)','(h)'};
labs3 = {'(i)','(j)','(k)','(l)'};

% Loop over each figure
for m=1:length(NPITHRESHOLD_ON0)
    NPIthreshold_on = NPITHRESHOLD_ON0(m);
    NPIthreshold_off = NPITHRESHOLD_OFF0(m);
    for a=1:length(ALPHA100)
        alpha1 = ALPHA100(a);
        alpha2 = alpha1;
        beta1 = R0*(alpha1+gamma)/N;
        DMAX = 1E-20;
        
        if(m==1)
            figure(a + 10)
        else
            figure(a + 13)
        end
        
        % Setup figure
        clf
        set(gcf,'color','w')
        set(gcf,'PaperUnits','centimeters')
        xSize = 10*4/3; ySize = 10;
        xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        set(gcf,'Position',[10 100 xSize*50 ySize*50])
        
        % Loop over variation in the relative transmission rate
        for s0=1:length(BETA2MULT)
            beta2mult = BETA2MULT(s0);
            beta2 = beta1*beta2mult;
            
            filename = strcat('Data/variantevo_summary_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat');
            if(exist(filename,'file'))
                load(filename)
                
                DMAX = max(DMAX,max(median_deaths(:)));
                
                % Plot data
                subplot(3,length(BETA2MULT),s0)
                imagesc(median_variant_recovered)
                set(gca,'ydir','normal')
                set(gca,'xtick',linspace(1,length(C),6),'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
                set(gca,'ytick',linspace(1,length(C),6),'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
                set(gca,'clim',[0,1])
                set(gca,'fontsize',10)
                temp = get(gca,'position');
                temp(1) = temp(1) - 0.05;
                set(gca,'position',temp);
                if(s0==length(BETA2MULT))
                    c=colorbar;
                    temp = get(c,'position');
                    temp(1) = temp(1) + 0.07;
                    set(c,'position',temp);
                    ylabel(c,{'Proportion infected','by variant'},'interpreter','latex','fontsize',14)
                end
                text(0,length(C)*1.1,labs1{s0},'interpreter','latex','fontsize',14)
                title(strcat('$\frac{\beta_v}{\beta_w}=',num2str(beta2mult),'$'),'interpreter','latex','fontsize',16)
                
                ax(s0) = subplot(3,length(BETA2MULT),s0+length(BETA2MULT));
                imagesc(median_deaths)
                set(gca,'ydir','normal')
                set(gca,'xtick',linspace(1,length(C),6),'xticklabel',{'0','0.2','0.4','0.6','0.8','1'})
                set(gca,'ytick',linspace(1,length(C),6),'yticklabel',{'0','0.2','0.4','0.6','0.8','1'})
                set(gca,'fontsize',10)
                temp = get(gca,'position');
                temp(1) = temp(1) - 0.05;
                set(gca,'position',temp);
                if(s0==1)
                    ylabel('Strength of cross immunity, $c$','interpreter','latex','fontsize',16)
                end
                if(s0==length(BETA2MULT))
                    c=colorbar;
                    temp = get(c,'position');
                    temp(1) = temp(1) + 0.07;
                    set(c,'position',temp);
                    ylabel(c,{'Total deaths','(per 100k)'},'interpreter','latex','fontsize',14')
                end
                text(0,length(C)*1.1,labs2{s0},'interpreter','latex','fontsize',14)
                
                subplot(3,length(BETA2MULT),s0+2*length(BETA2MULT))
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
                    x1=xlabel('Strength of NPIs, $r$','interpreter','latex','fontsize',16);
                    temp=get(x1,'position');
                    temp(1) = temp(1) + length(variant_emergence)*0.7;
                    set(x1,'position',temp);
                end
                if(s0==length(BETA2MULT))
                    c=colorbar;
                    temp = get(c,'position');
                    temp(1) = temp(1) + 0.07;
                    set(c,'position',temp);
                    ylabel(c,{'Variant emergence','probability'},'interpreter','latex','fontsize',14)
                end
                text(0,length(C)*1.1,labs3{s0},'interpreter','latex','fontsize',14)
            end
        end
        
        DMAX = 100*ceil(DMAX/100);
        for s0=1:length(BETA2MULT)
            axes(ax(s0))
            set(gca,'clim',[0,DMAX])
        end
        
%         % Save figures
%         if(m==1)
%             if(a==1)
%                 save2pdf('fig_S1.pdf');
%             else
%                 save2pdf('fig_S2.pdf');
%             end
%         else
%             if(a==1)
%                 save2pdf('fig_S4.pdf');
%             else
%                 save2pdf('fig_S5.pdf');
%             end
%         end
    end
end