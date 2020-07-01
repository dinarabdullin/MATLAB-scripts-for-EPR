function plotSpectrum(freq, spc, freqLim, spcLim, relFreq, saveFigure, filename)

    fig = figure('position', [0 0 600 600]);
    set(gcf, 'color', 'w');  
    
    spc = spc / max(spc);
    plot(freq, spc, 'Color', 'black', 'Linewidth', 1.5);
    set(gca, 'XLim', freqLim(1:2))
    set(gca, 'YLim', spcLim(1:2))
    set(gca, 'XTick', freqLim(1):freqLim(3):freqLim(2))
    set(gca, 'YTick', spcLim(1):spcLim(3):spcLim(2))
    if (relFreq)
        xlabel('{\it\nu}_{dd} (\nu_{0})');
    else
        xlabel('{\it\nu}_{dd}( MHz)');
    end
    ylabel('Spectrum (a.u.)');
    set(gca, 'fontsize', 16, 'TickDir', 'out');
    
	if (saveFigure)
        set(fig, 'PaperPositionMode', 'auto');
        saveas(fig, filename);
    end

end

