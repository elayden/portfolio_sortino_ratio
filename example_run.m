script_path = 'portfolio_sortino_ratio'; % insert path to folder 'portfolio_sortino_ratio'
addpath(script_path)

csv_dir = fullfile(script_path,'sample_portfolio');

%  [Sortino ratio, Sharpe ratio, total return, Downside Risk, SD, Drawdown]
outcomeWeights = [1,0,1,.09,0,1.5];

nRandom = 10000;  
minWeight = .05; maxWeight = 1; limitTickers = 10;
plotScatter = true; plotCorrs = false;
useParallel = true; 

% Run optimization:
results = portfolio_sortino_ratio(csv_dir, 'nRandom', nRandom, ...
    'outcomeWeights', outcomeWeights, 'minWeight', minWeight, ...
    'maxWeight', maxWeight,'limitTickers', limitTickers, ...
    'plotScatter', true, 'plotCorrs', plotCorrs,'parallel',useParallel) %#ok

% View optimized weights:
results.weights

% View average stats across years:
results.stats

% View returns for optimized portfolio compared to individual assets across
% years:
results.returns