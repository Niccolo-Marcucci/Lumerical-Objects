function nicePlot()
    set(gca,'FontSize',16,'LineWidth',1.5,'FontName','Arial',...
            'XGrid','on','YGrid','on')
    
    curves=findobj('type','line');
    for i=1:length(curves)
        set(curves(i),'linewidth',2)
    end
end
    