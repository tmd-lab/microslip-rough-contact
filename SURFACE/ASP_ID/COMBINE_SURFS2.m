function results = COMBINE_SURFS2(top, bot, combineSettings)
% Function for combining processed details from Two different surfaces into
% properties to be used in the dynamic simulations.
%
% Inputs:
%   top - surface processed for the top
%   bot - surface processed for the bottom
%   combineSettings - structure with all settings
%           .IQR_Factor   - multiple of IQR used to define outliers in the probability
%                distribution. 
%           .drawPlots    - true false to plot some different things
%           .NumInterpPoints - number of points to save probability at
%                                   distribution at.
%           .maxGaps      - maximum number of gaps to fit the probability
%                           with. Using all points can fill up memory and
%                           cause it to die.
%           .maxHistogram - maximum number of points for histogram plots if
%                           drawPlots = true
%
% Outputs:
%   results - structure with relevant fields:
%       Rprime - large principal radius of combined surfaces
%       Rpprime - small principal radius of combined surfaces
%       area_density - average area density of asperities on the two
%                       surfaces
%       zmax - maximum gap fit in the probability distribution, minimum is
%               zero.
%       normzinterp - values of z that can be used to interpolate probability
%                   distrubtion. - range of [0, 1]
%       pzinterp - values of probability distribution at zinterp
%
% NOTES:
%   1. Assumes that both surfaces are machined the same way such that the
%   statistics are shared. Mesoscale is assumed to be eliminated before
%   this function is called.
    
    %% Options (some of them)
    
    IQR_Factor = combineSettings.IQR_Factor;
    NumInterpPoints = combineSettings.NumInterpPoints;
    drawPlots = combineSettings.drawPlots;

    %% Straight Forward Asperity Stats
    
    % Reduce both data sets to only the relevant points
    mask_top = top.valid;
    mask_bot = bot.valid;
    
    median_r = median([top.radii(:, mask_top), bot.radii(:, mask_bot)], 2);
    
    % Asumes that two surfaces are statistically identical
    results.Rprime = median_r(2)/2;  %(1/median_r(2) + 1/median_r(2))^(-1);
    results.Rpprime = median_r(1)/2; %(1/median_r(1) + 1/median_r(1))^(-1);
    
    results.area_density = (sum(mask_top) + sum(mask_bot)) / (top.area + bot.area);
                            
    %% Gather Asperities
    
    top_tmp.z = top.xyzloc(mask_top, 3);
    bot_tmp.z = bot.xyzloc(mask_bot, 3);
    
    % Data structure for fitting plane and heights
    dat = {top_tmp, bot_tmp};
    
    
    %% Gap Distribution
    
    % Match the median height value between the two interfaces to combine
    % into one height distribution
    all_heights = [dat{1}.z - median(dat{1}.z); 
                    dat{2}.z - median(dat{2}.z)];
    
    % all of both interfaces to all of both interfaces
    gaps = -all_heights - all_heights';


    % eliminate outliers
    quants_z = quantile(gaps(:), 3);
    
    z_min = quants_z(1) - IQR_Factor*(quants_z(3) - quants_z(1));
    z_max = quants_z(3) + IQR_Factor*(quants_z(3) - quants_z(1));
    
    %save full gap data, with zero at the minimum used for the pdf
%     results.allGaps = gaps(:)-z_min;

    keep_asp = gaps >= z_min & gaps <= z_max;
    
    gaps = gaps(keep_asp);

    norm_gap = (gaps-z_min)/(z_max - z_min);
    
    if(length(norm_gap) > combineSettings.maxGaps)
        
       % Sort and then reduce
       norm_gap = sort(norm_gap);
       
       norm_gap = norm_gap(round(linspace(1, length(norm_gap), combineSettings.maxGaps)));
        
    end
    
    % Empirical fit of distrubtion + setting bounds
    emp_fit_flat = fitdist(norm_gap,'Kernel','Kernel','epanechnikov', 'Support', [0, 1]);
    
    % Other data to save
    results.z_max = z_max - z_min;
    
    results.normzinterp = linspace(0, 1, NumInterpPoints);
    results.pzinterp = pdf(emp_fit_flat,results.normzinterp);
    
    
    
    %% Plotting Options
    if(drawPlots)
                
        % Plot Gap Histogram and fit it
        set(groot, 'defaultAxesTickLabelInterpreter','default');  %Tex
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaultTextInterpreter','latex');

        
        %%%%%%%% Nice Plot for Paper
        unit_convert = 1e6;
        figure;
        hold on;
        h = histogram((gaps-z_min)*unit_convert);
        
        h.NumBins = 75;
        
        ylabel('Count');
        yrange = ylim;
        
        yyaxis right;
        
        xplot = linspace(0, 1, 100);
        yplot_emp = pdf(emp_fit_flat,xplot);
        plot(unit_convert*xplot*(z_max - z_min), yplot_emp, '-', 'LineWidth', 3);

        ylabel('Empirical PDF Value');

        ylim(yrange/length(gaps)*length(h.BinCounts))
        
        xlabel('Gap [$\mu$m]');
        
        set(gca, 'fontsize', 14);
        set(gcf, 'Renderer', 'painters');
        drawnow;
        
%         print('../PLOTS/OUTPUTS/paper_ltw_distribution.eps', '-depsc', '-r600');
        title('Reduced Gaps');
        
        
        %%%%%%%%%%%% Alpha Plot for Paper
        figure;
        hold on;
        alpha = [top.alpha(mask_top), bot.alpha(mask_bot)];
        alpha1 = alpha + pi*(alpha < 0);
        edges = linspace(0, pi, 32);
        
        h = histogram(alpha1, edges);
        xlim([0, pi]);
                
        ylabel('Count');
        xlabel('$\alpha$ [rad]')
        
        set(gca, 'fontsize', 14);
        set(gcf, 'Renderer', 'painters');
        drawnow;
        
%         print('../PLOTS/OUTPUTS/paper_alpha_distribution.eps', '-depsc', '-r600');
%         keyboard
        
        
        %%%%%%%%%%%% Alpha Plot for Paper
        figure;
        hold on;
        alpha = [top.alpha(mask_top), bot.alpha(mask_bot)];
        edges = linspace(-pi/2, pi/2, 32);
        
        h = histogram(alpha, edges);
        xlim([-pi/2, pi/2]);
                
        ylabel('Count');
        xlabel('$\alpha$ [rad]')
        
        set(gca, 'fontsize', 14);
        set(gcf, 'Renderer', 'painters');
        drawnow;
        
%         print('../PLOTS/OUTPUTS/paper_alpha_distribution.eps', '-depsc', '-r600');
%         print('../PLOTS/OUTPUTS/paper_alpha_distribution.svg', '-dsvg', '-r600');
%         addpath('../../export_fig/');
%         set(gcf, 'Color', 'w')
%         export_fig('../PLOTS/OUTPUTS/paper_alpha_distribution.png', '-dpng', '-r600'); %Needed to keep the ticks on the y axis, need png so no lines through colors.
%         keyboard
        
        
    end
    
end