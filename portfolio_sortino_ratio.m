function [results] = portfolio_sortino_ratio(csv_dir, varargin)
% RESULTS = PORTFOLIO_SORTINO_RATIO(CSV_DIR, varargin)
% 
% *Info:
%   This function optimizes portfolio weights based on a user-specified 
%   weighted linear combination of the Sortino ratio, Sharpe ratio, average 
%   total return, average downside risk, average standard deviation of 
%   returns, and max drawdown. The basic idea is to provide a directory as 
%   input 'csv_dir'. This folder should contain .csv files of historical 
%   data for each ticker desired to comprise part of a portfolio. The 
%   function will then return an optimized weighting scheme based on user's
%   criteria. It will also output historical performance data for the 
%   portfolio, an alternative equally-weighted portfolio, and for each 
%   individual ticker. Furthermore, you can choose to plot a matrix showing 
%   correlations among the individual assets, as well as a 3D point cloud 
%   showing where the optimized portfolio falls among other randomly 
%   generated portfolios on the highest-weighted dimensions. The function 
%   was designed using data from Yahoo Finance (https://finance.yahoo.com/) 
%   but should work with other data sources provided the formatting is 
%   similar. Try optimizing the sample portfolio included (with data from 
%   2006-2018) to gain a better understanding of use cases. (See example
%   below.) 
% 
% *Note: you must keep ^IRX.csv (13 week treasury bill rates) in the util 
%   folder up-to-date for risk free rate data. Download historical data at
%   <https://finance.yahoo.com/quote/%5EIRX/history?>. If other data in
%   your portfolio is more recent than ^IRX.csv, that additional data will
%   be discarded automatically. Also, statistics will only go back as far
%   as the most recent oldest data across all tickers.
% 
% *Cautionary note: your results may be skewed if your data does not go 
%   back sufficiently far and includes only one portion of a market cycle 
%   (e.g., all bull market, no recessions). To gain insight into 
%   performance over the whole market cycle, try to include data going back 
%   to 2008 or 2000, if not longer. When this is not possible, note that 
%   assets with high volatility and high annual returns during a bull 
%   market will often be the same assets that sustain the largest losses 
%   during major market corrections. In the case of major recessions these 
%   losses can sometimes exceed -40% in a calendar year.
% 
% *Disclaimer: This open-source research tool is not intended to provide 
%   investment advice. It is intended only for informational purposes, and 
%   the user is not recommended to use the tool to make actual investment 
%   decisions. Seek a duly licensed professional for investment advice.
% 
% Inputs:
% 
%   'csv_dir',        a full path which contains .csv files of financial 
%                     data; files should be titled using the ETF/stock 
%                     ticker name
% 
% Optional Name-Value Pair Arguments:
% 
%   'nRandom',        "optimization" uses the simple method of iterating 
%                     over a large number of random portfolio weights and
%                     then selecting the best based on given criteria. The
%                     larger this number, the better the selection that
%                     will likely be made. (default:  10,000)
% 
%   'outcomeWeights', a 1x6 vector containing weights which
%                     correspond to optimization preferences for
%                     [Sortino ratio, Sharpe ratio, total return, 
%                     Downside Risk, SD, Drawdown] (i.e., setting this 
%                     parameter allows one to gear the portfolio 
%                     optimization toward a more aggressive or conservative 
%                     portfolio strategy; weights need not sum to 1 
%                     (default: [1,0,1,0,0,1])
% 
%   'minWeight',      numeric value (range: 0-1) specifying the minimum 
%                     weight that any ticker may receive (default: 0)
% 
%   'maxWeight',      numeric value (range: 0-1) specifying the maximum 
%                     weight that any ticker may receive (default: 1)
% 
%   'limitTickers',   numeric value (range: 0 to nTickers) specifying the
%                     number of tickers to include in the optimized
%                     portfolio (i.e., can be a subset of the max number
%                     possible) (default:  include all tickers)
% 
%   'plotScatter',    boolean denoting whether to plot a 3D scatter plot
%                     showing the optimized portfolio among
%                     randomly permuted portfolios in terms of the most 
%                     highly weighted dimensions (default: true)
% 
%   'plotCorrs',      boolean denoting whether to plot an imagesc plot
%                     showing the correlations between the returns of each
%                     ticker across all years (default: true)
% 
%   'parallel',       boolean denoting whether to use parallel processing
%                     (default: false)
% 
% Outputs:
% 
%   'results',        a structure containing:
% 
%       'results.weights'
%                     a table showing tickers with optimized weights
% 
%       'results.stats'
%                     a table of average Sortino, Sharpe, total return, 
%                     downsideRisk, SD, and max drawdown across years for 
%                     the optimized portfolio, an equal weights portfolio, 
%                     and ranges across permutations for each measure
% 
%       'results.sortino'
%                     a table displaying Sortino ratio by year for both 
%                     the optimized portfolio, an equally weighted 
%                     portfolio, and for each individual ticker; tickers
%                     are sorted in descending order of average Sortino
% 
%       'results.sharpe'
%                     a table displaying Sharpe ratio by year for both 
%                     the optimized portfolio, an equally weighted 
%                     portfolio, and for each individual ticker; tickers
%                     are sorted in descending order of average Sharpe
% 
%       'results.returns'
%                     the same as (.sharpe) but for total return; tickers
%                     are sorted in descending order of average returns
% 
%       'results.downsideRisk'
%                     the same as (.sharpe) but for downside risk; tickers 
%                     are sorted in ascending order of average downside
%                     risk
% 
%       'results.SD'
%                     the same as (.sharpe) but for SD; tickers are 
%                     sorted in ascending order of average SD
% 
% Example function call:
% 
%   results = portfolio_sharpe_ratio(folderPath,...
%       'nRandom',10000,'outcomeWeights',[.3,0,.3,0,0,.4],...
%       'minWeight',.01,'maxWeight',.9,'limitTickers',10,'plotScatter',1,...
%       'plotCorrs',0,'parallel',true);
% 
%   -this will use 10,000 randomly permuted portfolio weights
%   -weights minimizing drawdown as a more important evaluation criterion 
%      than maximizing Sortino ratio or total return (both are weighted .3) 
%   -will not weight any portfolio asset less than 1%
%   -will not weight any portfolio asset more than 90%
%   -will limit the number of assets in the optimized portfolio to 10
%   -will plot the 3D scatter plot but not the corrlation matrix
%   -will use parallel processing
% 
% Author:  Elliot Layden, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hWait1 = waitbar(.5,'Please wait...');

% Identify Function Path and Add Helper Scripts:
script_fullpath = mfilename('fullpath'); 
[script_path,~,~] = fileparts(script_fullpath); 
addpath(genpath(script_path))

% Get Inputs:
inputs = varargin;
parsed = struct('nRandom', 10000, 'outcomeWeights', [1,0,1,0,0,1], ...
    'minWeight', 0, 'maxWeight', 1,'limitTickers', 0, 'plotScatter', true,...
    'plotCorrs', true,'parallel',false);
% Column of internal cells == OR, Row of internal cells == AND
input_types = {{'numeric'},{'vector'},{'numeric'},{'numeric'},{'numeric'},...
    {'logical';'numeric'},{'logical';'numeric'},{'logical';'numeric'}};  
parsed = getInputs(inputs, parsed, input_types);
nRandom = parsed.nRandom; outcomeWeights = parsed.outcomeWeights;
limitTickers = parsed.limitTickers;
minWeight = parsed.minWeight; maxWeight = parsed.maxWeight;
plotScatter = parsed.plotScatter; plotCorrs = parsed.plotCorrs;
useParallel = parsed.parallel;
if nRandom <=0; nRandom = 10000; end
if length(outcomeWeights)~=6
    warning('''outcomeWeights'' should be a 1x6 vector. Please check your inputs.')
    outcomeWeights = [1,0,1,0,0,1]; 
end
outcomeWeights = outcomeWeights./sum(outcomeWeights);

% Initialize Results:
results = struct('weights',[],'stats',[],'sortino',[],'sharpe',[],...
    'returns',[],'downsideRisk',[],'SD',[]);

% Extract portfolio data:
warning('off','all'); data = extract_portfolio(csv_dir); warning('on','all')

% Get risk-free rates (10-year treasury rates):
warning('off','all')
riskFreeData = readtable(fullfile(script_path,'util','^IRX.csv'));
warning('on','all')
riskFreeData.Date = datenum(riskFreeData.Date);
[~,ixUse,ixRiskFree] = intersect(data(1).date,riskFreeData.Date);
riskFreeData = riskFreeData(ixRiskFree,:);
annualized = str2double(riskFreeData.Close);
riskFree = (( 1 + annualized/100 ).^(1/252)-1); % data starts as whole-number percentages
ixDelete = isnan(riskFree); riskFree(ixDelete) = [];
dates = year(riskFreeData.Date); dates(ixDelete) = []; 
years = unique(dates); nYears = length(years);

% Limit portfolio data to what's available for riskFree data:
for ix = 1:length(data)
    data(ix).date = data(ix).date(ixUse);
    data(ix).raw = data(ix).raw(ixUse);
    data(ix).return = data(ix).return(ixUse);
    data(ix).volume = data(ix).volume(ixUse);
    data(ix).date(ixDelete) = [];
    data(ix).raw(ixDelete) = [];
    data(ix).return(ixDelete) = [];
    data(ix).volume(ixDelete) = [];
end

% Extract all tickers:
nRows = size(data(1).return,1); nTickers = length(data);
returns = zeros(nRows,nTickers); raw = returns; tickers = cell(1, nTickers);
for ix = 1:nTickers
   returns(:,ix) = data(ix).return; 
   raw(:,ix) = data(ix).raw; 
   tickers{ix} = data(ix).ticker;
end

% Input Checking:
if limitTickers >= nTickers; limitTickers = 0; end
if minWeight < 0; minWeight = 0; end
if maxWeight > 1; maxWeight = 1; end

axesLabels = {'Sortino','Sharpe','Return','DownsideRisk','SD','Drawdown'};

% Permuted weights:
delete(hWait1); 
if useParallel
    hCheck = gcp('nocreate');
    if isempty(hCheck) || ~hCheck.Connected
        hPool = parpool; 
    end
else
    hWait = waitbar(0,'Processing...');
end
if limitTickers > 0
    weights = randfixedsum(limitTickers,nRandom,1,minWeight,maxWeight)';
    useTickers = zeros(nRandom, limitTickers);
else
    weights = randfixedsum(nTickers,nRandom,1,minWeight,maxWeight)';
    useTickers = repmat((1:nTickers), nRandom, 1);
end

meanSharpe = zeros(nRandom,1); meanSortino = zeros(nRandom,1);
meanReturn = zeros(nRandom,1); meanSD = zeros(nRandom,1);
maxDrawdown = zeros(nRandom,1); meanDownsideRisk = zeros(nRandom,1);
if useParallel
    fprintf('Parallel computing progress:\n');
    fprintf(['\n' repmat('.',1,nRandom/1000) '\n\n']);
end
if useParallel
    parfor perm = 1:nRandom
        if rem(perm,1000)==0
            fprintf('\b|\n'); % prints a tick at start of each parallel iteration
        end
        if limitTickers > 0
            useTickers(perm,:) = randperm(nTickers, limitTickers);
        end
        sharpe = zeros(1,nYears); sortino = sharpe; yearlyReturns = sharpe; SD = sharpe; downsideRisk = sharpe;
        for ix = 1:nYears
            [sharpe(ix), sortino(ix), yearlyReturns(ix), SD(ix), downsideRisk(ix)] = calculate_metrics(...
                returns(dates==years(ix),useTickers(perm,:)), ...
                weights(perm,:), riskFree(dates==years(ix))); %#ok
        end
        meanSortino(perm) = geomean(sortino +10);
        meanSortino(perm) = meanSortino(perm) -10;
        meanSharpe(perm) = geomean(sharpe +10); 
        meanSharpe(perm) = meanSharpe(perm) -10;
        meanReturn(perm) = geomean(yearlyReturns + 10);
        meanReturn(perm) = meanReturn(perm) -10;
        meanSD(perm) = geomean(SD);
        meanDownsideRisk(perm) = geomean(downsideRisk);
        maxDrawdown(perm) = min(yearlyReturns);
    end
else
    for perm = 1:nRandom
        if rem(perm,1000)==0
            waitbar(perm/nRandom, hWait);
        end
        if limitTickers > 0
            useTickers(perm,:) = randperm(nTickers, limitTickers);
        end
        sharpe = zeros(1,nYears); sortino = sharpe; yearlyReturns = sharpe; SD = sharpe; downsideRisk = sharpe;
        for ix = 1:nYears
            [sharpe(ix), sortino(ix), yearlyReturns(ix), SD(ix), downsideRisk(ix)] = calculate_metrics(...
                returns(dates==years(ix),useTickers(perm,:)), ...
                weights(perm,:), riskFree(dates==years(ix)));
        end
        meanSortino(perm) = geomean(sortino +10);
        meanSortino(perm) = meanSortino(perm) -10;
        meanSharpe(perm) = geomean(sharpe +10); 
        meanSharpe(perm) = meanSharpe(perm) -10;
        meanReturn(perm) = geomean(yearlyReturns + 10);
        meanReturn(perm) = meanReturn(perm) -10;
        meanSD(perm) = geomean(SD);
        meanDownsideRisk(perm) = geomean(downsideRisk);
        maxDrawdown(perm) = min(yearlyReturns);
    end
end
if useParallel && (isempty(hCheck) || ~hCheck.Connected)
    delete(hPool)
else
    delete(hWait)
end

% Calculate Mahalanobis Distance of each permuted portfolio to ideal:
permData = [meanSortino, meanSharpe, meanReturn, meanDownsideRisk, meanSD, maxDrawdown];
permDataZ = zscore(permData);
permDataZ = permDataZ - repmat(min(permDataZ),nRandom,1);
permDataZ = permDataZ ./ repmat(max(permDataZ),nRandom,1);
% best = repmat([max(permDataZ(:,1)), min(permDataZ(:,2)), max(permDataZ(:,3))],nRandom,1);
% dist = sqrt(sum(repmat(outcomeWeights.*[1,-1,1],nRandom,1).*(permDataZ - repmat(best,nRandom,1)).^2,2));
dist = sqrt(sum(repmat(outcomeWeights,nRandom,1).*(permDataZ - repmat([1,1,1,0,0,1],nRandom,1)).^2,2));
[~, ix_best] = min(dist); % find min. distance

compute_results(ix_best);

function compute_results(ix_best)
    
    bestWeights = weights(ix_best,:); % find best weights
    [bestWeights_sort, ixWeights] = sort(bestWeights,'descend');
    
    % Add to results struct:
    results.weights = cell2table(num2cell(round(bestWeights_sort*100,3)));
    tickers1 = tickers(useTickers(ix_best,:)); tickers1 = tickers1(ixWeights);
    results.weights.Properties.VariableNames = tickers1;
    results.weights.Properties.RowNames = {'% weight:'};

    % Calculate yearly for optimized weights:
    sharpe_opt = zeros(1,nYears); sortino_opt = sharpe_opt;
    yearlyReturns_opt = sharpe_opt; SD_opt = sharpe_opt; downsideRisk_opt = sharpe_opt;
    for i = 1:nYears
        [sharpe_opt(i), sortino_opt(i), yearlyReturns_opt(i), SD_opt(i), downsideRisk_opt(i)] = calculate_metrics(...
            returns(dates==years(i),useTickers(ix_best,:)), bestWeights, riskFree(dates==years(i)));
    end    
    % Calculate for equal weighting:
    equalWeights = ones(1,nTickers)./nTickers; 
    sharpe_eq = zeros(1,nYears); sortino_eq = sharpe_eq;
    yearlyReturns_eq = sharpe_eq; SD_eq = sharpe_eq; downsideRisk_eq = sharpe_eq;
    for i = 1:nYears
        [sharpe_eq(i), sortino_eq(i), yearlyReturns_eq(i), SD_eq(i), downsideRisk_eq(i)] = calculate_metrics(...
            returns(dates==years(i),:), equalWeights, riskFree(dates==years(i)));
    end
    % Calculate stats for individual tickers:
    sharpe_ind = zeros(nTickers, nYears); sortino_ind = sharpe_ind;
    yearlyReturns_ind = sharpe_ind; SD_ind = sharpe_ind; downsideRisk_ind = sharpe_ind;
    for i = 1:nTickers
        for j = 1:nYears
            [sharpe_ind(i,j), sortino_ind(i,j), yearlyReturns_ind(i,j), SD_ind(i,j), downsideRisk_ind(i,j)] = calculate_metrics(...
                returns(dates==years(j),i), 1, riskFree(dates==years(j)));
        end
    end

    % Save table output 'stats':
    results.stats = cell2table(num2cell(round([...
        geomean(sortino_opt + 10), geomean(sharpe_opt +10), ...
        geomean(yearlyReturns_opt + 10), geomean(downsideRisk_opt), ...
        geomean(SD_opt) min(yearlyReturns_opt);...
        geomean(sortino_eq + 10), geomean(sharpe_eq +10), ...
        geomean(yearlyReturns_eq +10), geomean(downsideRisk_eq), ...
        geomean(SD_eq), min(yearlyReturns_eq);...
        max(meanSortino), max(meanSharpe), max(meanReturn), ...
        min(meanDownsideRisk), min(meanSD), max(maxDrawdown);...
        min(meanSortino), min(meanSharpe), min(meanReturn), ...
        max(meanDownsideRisk), max(meanSD), min(maxDrawdown)],3)));
    results.stats.Var1(1:2) = results.stats.Var1(1:2) - 10;
    results.stats.Var2(1:2) = results.stats.Var2(1:2) - 10;
    results.stats.Var3(1:2) = results.stats.Var3(1:2) - 10;
    results.stats.Properties.VariableNames = axesLabels;
    results.stats.Properties.RowNames = {'Optimized','EqualWeights','Best','Worst'};

    yearLabels = sprintfc('%d',years);
    for i = 1:length(yearLabels)
        yearLabels{i} = ['y',yearLabels{i}];
    end
    rowLabels = [{'Optimized','EqualWeights'},tickers];

    % Save table output 'sortino':
    sortinoMat = [sortino_opt; sortino_eq; sortino_ind];
    [~, sortinoMat_ix] = sort(geomean(sortinoMat + 10,2),'descend');
    results.sortino = cell2table(num2cell(round(sortinoMat(sortinoMat_ix,:),3)));
    results.sortino.Properties.VariableNames = yearLabels;
    results.sortino.Properties.RowNames = rowLabels(sortinoMat_ix);

    % Save table output 'sharpe':
    sharpeMat = [sharpe_opt; sharpe_eq; sharpe_ind];
    [~, sharpeMat_ix] = sort(geomean(sharpeMat + 10,2),'descend');
    results.sharpe = cell2table(num2cell(round(sharpeMat(sharpeMat_ix,:),3)));
    results.sharpe.Properties.VariableNames = yearLabels;
    results.sharpe.Properties.RowNames = rowLabels(sharpeMat_ix);

    % Save table output 'returns':
    yearlyReturnsMat = [yearlyReturns_opt; yearlyReturns_eq; yearlyReturns_ind];
    [~, yearlyReturns_ix] = sort(geomean(yearlyReturnsMat + 10,2),'descend');
    results.returns = cell2table(num2cell(round(yearlyReturnsMat(yearlyReturns_ix,:),3)));
    results.returns.Properties.VariableNames = yearLabels;
    results.returns.Properties.RowNames = rowLabels(yearlyReturns_ix);

    % Save table output 'DownsideRisk':
    DownsideRiskMat = [downsideRisk_opt; downsideRisk_eq; downsideRisk_ind];
    [~, downsideRisk_ix] = sort(geomean(DownsideRiskMat,2),'ascend');
    results.downsideRisk = cell2table(num2cell(round(DownsideRiskMat(downsideRisk_ix,:),3)));
    results.downsideRisk.Properties.VariableNames = yearLabels;
    results.downsideRisk.Properties.RowNames = rowLabels(downsideRisk_ix);

    % Save table output 'SD':
    SDMat = [SD_opt; SD_eq; SD_ind];
    [~, SD_ix] = sort(geomean(SDMat,2),'ascend');
    results.SD = cell2table(num2cell(round(SDMat(SD_ix,:),3)));
    results.SD.Properties.VariableNames = yearLabels;
    results.SD.Properties.RowNames = rowLabels(SD_ix);
    
    assignin('base','results',results)
end

% Plot Sharpe x Return x SD
rotate_on = true; zoom_on = false; crosshair_on = false;
targetSize = 200;
if plotScatter
    if nRandom > 10000
        sparseIx = randperm(nRandom,10000);
    else
        sparseIx = 1:nRandom;
    end
    [~,ix_axes] = sort(outcomeWeights,'descend');
    fig_title = [axesLabels{ix_axes(1)},' x ', axesLabels{ix_axes(2)},' x ',...
        axesLabels{ix_axes(3)},' x ',axesLabels{ix_axes(4)},' x ',axesLabels{ix_axes(5)}];
    hfig = figure('Name',fig_title,'units','norm','Color',[1,1,1],'Position',...
        [.061,.116,.861,.806],'NumberTitle','off','MenuBar','none'); 
    tool_menu = uimenu(hfig,'Label','Tools');
        h_tools(1) = uimenu(tool_menu,'Label','Rotate','Checked','on','Callback',{@change_tool,1});
        h_tools(2) = uimenu(tool_menu,'Label','Zoom','Callback',{@change_tool,2});
        h_tools(3) = uimenu(tool_menu,'Label','Select New','Callback',{@change_tool,3});
    hScatter = scatter3(permData(sparseIx,ix_axes(1)),permData(sparseIx,ix_axes(2)),...
        permData(sparseIx,ix_axes(3)),20,'o','fill'); hold on;
    mainAx = gca; set(mainAx,'FontSize',14)
    % Set point colors to outcomeWeights #4:
    cMax = max(permData(sparseIx,ix_axes(4))); cMin = min(permData(sparseIx,ix_axes(4)));
    useColors = round((permData(sparseIx,ix_axes(4))-cMin) ./ (cMax - cMin) .* 359)+1;
    colorMat = jet(360); colormap(colorMat);
    hScatter.CData = colorMat(useColors,:); set(gca,'clim',[cMin,cMax])
    hColor = colorbar; title(hColor,axesLabels{ix_axes(4)},'FontSize',12,'FontWeight','normal');
    % Add optimized portfolio:
    hCurrPoint1 = scatter3(permData(ix_best,ix_axes(1)), ...
        permData(ix_best,ix_axes(2)), permData(ix_best,ix_axes(3)),...
        targetSize,'x','k');
    hCurrPoint2 = scatter3(permData(ix_best,ix_axes(1)), ...
        permData(ix_best,ix_axes(2)), permData(ix_best,ix_axes(3)),...
        targetSize,'o','k');
    xlabel(axesLabels{ix_axes(1)}); ylabel(axesLabels{ix_axes(2)}); zlabel(axesLabels{ix_axes(3)}); rotate3d
    % Set sizes to outcomeWeights #5:
    minSize = 5; maxSize = 60;
    sizes = permData(sparseIx,ix_axes(5)) - min(permData(sparseIx,ix_axes(5)));
    sizes = (sizes./range(sizes)).*(maxSize-minSize) + minSize;
    hScatter.SizeData = sizes;
    leg_ax = axes('Parent',hfig,'Visible','off'); hold(leg_ax,'on'); h_leg = zeros(1,2);
    for ix = 1:2; h_leg(ix) = scatter(1,1,1,'filled','MarkerFaceColor','b','Parent',leg_ax); end
    
%     [h_legend, icons] = legend(mainAx, h_leg, 'Location','northwest','Labels',...
%         {num2str(round(min(permData(sparseIx,ix_axes(5))),4)), ...
%         num2str(round(max(permData(sparseIx,ix_axes(5))),4))},'FontName','Arial','FontSize',14);

    axes(mainAx)
    [h_legend, icons] = legend(h_leg, {num2str(round(min(permData(sparseIx,ix_axes(5))),4)), ...
        num2str(round(max(permData(sparseIx,ix_axes(5))),4))}, 'Location','northwest',...
        'FontName','Arial','FontSize',14);
    
    icons(3).Children.MarkerSize = 2; icons(4).Children.MarkerSize = 8;
    axes(mainAx)
    set(get(h_legend,'title'),'String',axesLabels{ix_axes(5)},...
        'FontSize',12,'FontWeight','normal','Visible','on')
    h_legend.Title.NodeChildren.Position= [0.5 1.5 0];
    set(leg_ax,'Visible','off'); 
end

% Plot ticker correlation matrix:
if plotCorrs
    nTickers = length(data); r = corr(returns); r(eye(nTickers)==1) = 0;
    fig_title = 'Ticker Correlation Matrix';
    figure('Name',fig_title,'units','norm','Color',[1,1,1],'Position',...
        [.236,.106,.554,.806],'NumberTitle','off','MenuBar','none'); 
    imagesc(r); 
    set(gca,'XTick',1:nTickers,'XTickLabels',tickers,'YTick',1:nTickers,...
        'YTickLabels',tickers,'XTickLabelRotation',90,'FontSize',12)
    colormap('jet'); colorbar; 
end

%% Callbacks:
function change_tool(~, ~, which_tool)
     switch which_tool
        case 1 % Rotate
            if ~rotate_on
                zoom off; rotate3d on; 
                rotate_on = true; zoom_on = false; crosshair_on = false;
                for ixx = 1:3; h_tools(ixx).Checked = 'off'; end
                h_tools(1).Checked = 'on';
            else
                rotate3d off; rotate_on = false; h_tools(1).Checked = 'off';
            end
                
        case 2 % Zoom
            if ~zoom_on
                zoom on; rotate3d off; 
                zoom_on = true; rotate_on = false; crosshair_on = false;
                for ixx = 1:3; h_tools(ixx).Checked = 'off'; end
                h_tools(2).Checked = 'on';
            else
                zoom off; zoom_on = false; h_tools(2).Checked = 'off';
            end
         case 3 % Cross-hair
             if ~crosshair_on
                zoom off; rotate3d off; 
                zoom_on = false; rotate_on = false; crosshair_on = true;
                for ixx = 1:3; h_tools(ixx).Checked = 'off'; end
                h_tools(3).Checked = 'on';
                set(hfig,'WindowButtonDownFcn',@cursor_click_callback,'pointer','crosshair'); 
             else
                crosshair_on = false; h_tools(3).Checked = 'off';
                set(hfig,'WindowButtonDownFcn',[],'pointer','arrow');
             end
     end
end

function cursor_click_callback(~,~,~)
    if crosshair_on
        % Find clicked point (Inspired by Click3dPoint by Babak Taati
        % (https://www.mathworks.com/matlabcentral/fileexchange/7594-click3dpoint)):
        pt = get(mainAx,'CurrentPoint'); % Get mouse location
        camPos = get(mainAx, 'CameraPosition'); % camera position
        camTgt = get(mainAx, 'CameraTarget'); % where the camera is pointing to
        camDir = camPos - camTgt; % camera direction
        camUpVect = get(mainAx, 'CameraUpVector'); % camera 'up' vector
        % Build an orthonormal frame based on the viewing direction and the up vector (the "view frame")
        zAxis = camDir/norm(camDir);    
        upAxis = camUpVect/norm(camUpVect); 
        xAxis = cross(upAxis, zAxis);
        yAxis = cross(zAxis, xAxis);
        rot = [xAxis; yAxis; zAxis]; % view rotation 
        % the point cloud represented in the view frame
        rotatedPointCloud = rot * permData(:,ix_axes(1:3))'; 
        % the clicked point represented in the view frame
        rotatedPointFront = rot * pt' ;
        % find the nearest neighbour to the clicked point 
        ix_best = dsearchn(rotatedPointCloud(1:2,:)',rotatedPointFront(1:2));
        selectedPoint = round(permData(ix_best,ix_axes),4);

        % Recompute results & adjust plot for selectedPoint:
        compute_results(ix_best);
        delete(hCurrPoint1); delete(hCurrPoint2);
        hCurrPoint1 = scatter3(permData(ix_best,ix_axes(1)), ...
            permData(ix_best,ix_axes(2)), permData(ix_best,ix_axes(3)), ...
            targetSize,'x','k');
        hCurrPoint2 = scatter3(permData(ix_best,ix_axes(1)), ...
            permData(ix_best,ix_axes(2)), permData(ix_best,ix_axes(3)), ...
            targetSize,'o','k');
        disp(['Results successfully recalculated for selected portfolio:  ',...
            axesLabels{ix_axes(1)},': ',num2str(selectedPoint(1)),';  ',...
            axesLabels{ix_axes(2)},': ',num2str(selectedPoint(2)),';  ',...
            axesLabels{ix_axes(3)},': ',num2str(selectedPoint(3)),';  ',...
            axesLabels{ix_axes(4)},': ',num2str(selectedPoint(4)),';  ',...
            axesLabels{ix_axes(5)},': ',num2str(selectedPoint(5))])
    end
end

end