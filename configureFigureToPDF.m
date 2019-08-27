function fig = configureFigureToPDF(fig,figName)
if exist('figName','var')
    set(fig,'RendererMode','manual','Renderer','painters',...
        'PaperOrientation','landscape','Color',[1,1,1],'Name',figName);
else
    set(fig,'RendererMode','manual','Renderer','painters',...
        'PaperOrientation','landscape','Color',[1,1,1]);
end
end