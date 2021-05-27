function fig_3

% fig_3.m
%
% Generates figure 3

% Variables
R = [0,0.5,0.75];

% Setup figure
figure(3)
clf
set(gcf,'color','w')
set(gcf,'PaperUnits','centimeters')
xSize = 14; ySize = 8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[10 100 xSize*50 ySize*50])

% Figure labels
labs1 = {'(a)','(b)','(c)'};
labs2 = {'(d)','(e)','(f)'};

% Load and plot data
for i=1:3
    r = R(i);
    
    load(strcat('Data/fig_3_data_',num2str(r),'.mat'))
    
    subplot(2,3,i)
    hold on
    for s=1:simtotal
        plot(0:MaxTime,(I1_OUT(:,s)+I21_OUT(:,s)+C_OUT(:,s))/N,'color',[0,0,1,0.25],'linewidth',1)
        plot(0:MaxTime,(I2_OUT(:,s)+I12_OUT(:,s)+C_OUT(:,s))/N,'color',[1,0,0,0.25],'linewidth',1)
    end
    if(i==1)
        legend('Wildtype','Variant')
    end
    set(gca,'fontsize',10)
    set(gca,'yscale','log')
    xlim([0,300])
    ylim([0.9/N,1])
    box on
    if(i==1)
        ylabel('Proportion infected','interpreter','latex','fontsize',16)
    end
    title(strcat('$r=',num2str(r),'$'),'interpreter','latex','fontsize',16)
    text(0,3,labs1{i},'interpreter','latex','fontsize',14)
    
    subplot(2,3,3+i)
    hold on
    for s=1:simtotal
        plot(0:MaxTime,(I2_OUT(:,s)+I12_OUT(:,s)+C_OUT(:,s))./(I1_OUT(:,s)+I2_OUT(:,s)+I12_OUT(:,s)+I21_OUT(:,s)+C_OUT(:,s)),'color',[0,0,0,0.35],'linewidth',1)
    end
    set(gca,'fontsize',10)
    set(gca,'yscale','log')
    xlim([0,300])
    ylim([0.9/N,1])
    box on
    set(gca,'ytick',logspace(log10(1/N),0,6))
    if(i==2)
        xlabel('Time (days)','interpreter','latex','fontsize',16)
    end
    if(i==1)
        ylabel('Variant frequency','interpreter','latex','fontsize',16)
    end
    text(0,3,labs2{i},'interpreter','latex','fontsize',14)
end

% Save figure
% save2pdf('fig_3.pdf')