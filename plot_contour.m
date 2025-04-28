function plot_contour(x, y, z, zlim, baseline, plotcolor, plottitle, laboff, bands, bandlabels, threshhold, testout)


    if isempty(testout)
        figure; colormap(plotcolor)
        contourf(x, y, z, 40, 'linecolor', 'none')
        set(gca,'clim',zlim); colorbar
        xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
        xline(0,'w','LineWidth',1)
        xline(baseline,'--w','LineWidth',1)
        yline(bands(3,1),'white','LineWidth',1); yline(bands(3,2),'white','LineWidth',1); % Theta
        yline(bands(4,1),'white','LineWidth',1); yline(bands(4,2),'white','LineWidth',1); % Beta
        t = text(x(1)+laboff,3,bandlabels(2)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,10,bandlabels(3)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,23,bandlabels(4)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,40,"low γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,65,"med γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,110,"high γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(plottitle,'FontSize',18);

    end

    if ~isempty(testout)
        figure; colormap(plotcolor); hold on
        contourf(x, y, z, 40, 'linecolor', 'none','FaceAlpha',.3)
        threshmap = normalize(testout,2) > threshhold; threshmap = z.*threshmap; threshmap(threshmap==0) = NaN;
        contourf(x, y, threshmap, 20, 'linecolor', 'none')
        set(gca,'clim',zlim); colorbar
        xlabel('Time (ms)'), ylabel('Frequencies (Hz)')
        xline(0,'w','LineWidth',1)
        xline(baseline,'--w','LineWidth',1)
        yline(bands(3,1),'white','LineWidth',1); yline(bands(3,2),'white','LineWidth',1); % Theta
        yline(bands(4,1),'white','LineWidth',1); yline(bands(4,2),'white','LineWidth',1); % Beta
        t = text(x(1)+laboff,3,bandlabels(2)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,10,bandlabels(3)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,23,bandlabels(4)); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,40,"low γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,65,"med γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        t = text(x(1)+laboff,110,"high γ"); t(1).Color = 'white'; t(1).FontSize = 14;
        xlabel('Time (s)'); ylabel('Frequency (Hz)'); title(plottitle,'FontSize',18);

    end

end