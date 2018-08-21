function [sharpe, sortino, totalReturn, portfolioSD, downsideRisk] = calculate_metrics(returns, weights, riskFreeReturn)

    useRiskAdjustedSortino = false; % if so, uses same risk adjustment as for Sharpe, rather than 0
    [nRows, nCols] = size(returns);
    
    % Geometric Sharpe Ratio: https://quant.stackexchange.com/questions/3607/should-i-use-an-arithmetic-or-a-geometric-calculation-for-the-sharpe-ratio
    % average exp(mean(log(1+returns))) - 1 is equiv. to geomean(returns)
    riskAdjustedReturns = returns - repmat(riskFreeReturn,1,nCols);
    compoundedAdjustedReturns = sum(repmat(weights,nRows,1).*log(1+riskAdjustedReturns),2);
    sharpe = mean(compoundedAdjustedReturns) ./ sum(weights.*std(log(1+riskAdjustedReturns))) .* sqrt(252);
    
    % Arithmetic:
%         SD = std(returns); sharpe = mean(riskAdjustedReturns) ./ SD .* sqrt(252); SD = SD .* sqrt(252); % annualized SD's

    % Expected Portfolio Return: https://quant.stackexchange.com/questions/25371/geometric-means-standard-deviation-and-sharpe-ratios
    logReturns = log(1+returns);
    totalReturn = (exp(mean(logReturns))).^nRows - 1;
    totalReturn = sum(totalReturn.*weights);
    
    % Portfolio SD:
    SD = std(logReturns).*sqrt(252); % annualized SD's
    wSD = weights .* SD;
    portfolioSD = sqrt(wSD * corr(logReturns) * wSD');
    
    % Sortino Ratio
    if useRiskAdjustedSortino
        downsideRisk = sum(weights .* sqrt( sum( min(0,log(1+riskAdjustedReturns)).^2 ) / nRows ));
        sortino = mean(compoundedAdjustedReturns) / downsideRisk * sqrt(252);
    else
        downsideRisk = sum(weights .* sqrt( sum( min(0,log(1+returns)).^2 ) / nRows ));
        sortino = mean(sum(repmat(weights,nRows,1).*log(1+returns),2)) / downsideRisk * sqrt(252);
    end
end