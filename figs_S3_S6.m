function figs_S3_S6

% figs_S3_S6.m
%
% Generates figures S3 and S6

% Fixed parameters
N = 1e5;
R0 = 3;
gamma = 1/5;
beta2mult = 2;

% Variables
ALPHA10 = gamma./[499,4];
NPITHRESHOLD_ON = [0,0.01];
NPITHRESHOLD_OFF = [0,0.002];
R = linspace(0,1,101);

% Figure labels
labs1 = {'(a)','(b)','(c)','(d)','(e)','(f)'};

% Loop over each figure

for k=1:2
    NPIthreshold_on = NPITHRESHOLD_ON(k);
    NPIthreshold_off = NPITHRESHOLD_OFF(k);
    if(k==1)
        figure(13)
    else
        figure(16)
    end
    
    % Setup figure
    clf
    set(gcf,'color','w')
    set(gcf,'PaperUnits','centimeters')
    xSize = 14; ySize = 8;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[10 100 xSize*50 ySize*50])
    
    for a0=1:length(ALPHA10)
        alpha1 = ALPHA10(a0);
        ALPHA2 = alpha1*[1,2,3];
        beta1 = R0*(alpha1+gamma)/N;
        beta2 = beta1*beta2mult;
    
        ymax = 0;
        
        % Loop over variation in the relative mortality rate
        for l=1:length(ALPHA2)
            alpha2 = ALPHA2(l);
            
            subplot(length(ALPHA10),length(ALPHA2),l + (a0-1)*3)
            hold on
            
            filename = strcat('Data/variantevo_mortality_',num2str(NPIthreshold_on),'_',num2str(NPIthreshold_off),'_',num2str(alpha1),'_',num2str(alpha2),'_',num2str(beta1),'_',num2str(beta2),'.mat');
            load(filename)
            
            % Plot data
            OUT = prctile(Dtotal(1,:,:),[25,75],3);
            patch([R fliplr(R)],[max(0,OUT(:,:,1)) fliplr(OUT(:,:,2))],'k','facealpha',0.3,'linestyle','none');
            plot(R,median(Dtotal(end,:,:),3),'color','k','linewidth',2)
            
            ymax = max(ymax,max(OUT(:,:,2)));
            
            box on
            if(a0==length(ALPHA10) && l==2)
                xlabel('Strength of NPIs, $r$','interpreter','latex','fontsize',16)
            end
            if(a0==length(ALPHA10) && l==1)
                y1=ylabel('Total deaths (per 100k)','interpreter','latex','fontsize',16);
                temp=get(y1,'position');
                temp(1) = temp(1) - 0.05;
                if(k==1)
                temp(2) = temp(2) + ymax;
                elseif(k==2)
                temp(2) = temp(2) + 1.2*ymax;
                end
                set(y1,'position',temp);
            end
            if(a0==1)
                if(alpha1==alpha2)
                    title(strcat('$\alpha_v=\alpha_w$'),'interpreter','latex','fontsize',16)
                else
                    title(strcat('$\alpha_v=',num2str(alpha2/alpha1),'\alpha_w$'),'interpreter','latex','fontsize',16)
                end
            end
        end
        
        ymax = 100*ceil(ymax/100);
        for l=1:length(ALPHA2)
            subplot(length(ALPHA10),length(ALPHA2),l + (a0-1)*3)
            ylim([0,ymax]);
            if(a0==length(ALPHA10))
                text(0,max(get(gca,'ylim'))*1.12,labs1{l + (a0-1)*3},'interpreter','latex','fontsize',14)
            else
                text(0,max(get(gca,'ylim'))*1.1,labs1{l + (a0-1)*3},'interpreter','latex','fontsize',14)
            end
        end
        
%         % Save figures
%         if(k==1)
%             save2pdf('fig_S3.pdf')
%         else
%             save2pdf('fig_S6.pdf')
%         end
    end
end