function [out] = extract_portfolio(csv_dir)
    % Find files:
    listing = dir(fullfile(csv_dir,'*.csv'));
    if isempty(listing)
        listing = dir(fullfile(csv_dir,'*.xlsx'));
    end
    % Extract tables:
    data = cell(1,length(listing));
    for i = 1:length(listing)
        dataTable = readtable(fullfile(csv_dir,listing(i).name));
        data{i} = dataTable;
    end
    % Get Ticker Names:
    tickers = cell(1,length(listing));
    for i = 1:length(listing)
        tickers{i} = listing(i).name(1:end-4);
    end
    % Check for presence of Close or AdjClose fields:
    for i = 1:length(data)
        colNames = data{i}.Properties.VariableNames;
        if sum(strcmpi(colNames,'AdjClose'))==0 && sum(strcmpi(colNames,'Close'))==0
            error(['ERROR:  Ticker ',tickers{i},' has no field ''AdjClose'' or ''Close''.'])
        end
    end
    % Check format of columnNames:
    for i = 1:length(data)
        colNames = data{i}.Properties.VariableNames;
        if sum(strcmpi(colNames,'Date'))
            data{i}.Properties.VariableNames{strcmpi(colNames,'Date')} = 'Date';
        end
        if strcmpi(colNames,'AdjClose')
            data{i}.Properties.VariableNames{strcmpi(colNames,'AdjClose')} = 'AdjClose';
        end
        if strcmpi(colNames,'Close')
            data{i}.Properties.VariableNames{strcmpi(colNames,'Close')} = 'Close';
        end
        if strcmpi(colNames,'volume')
            data{i}.Properties.VariableNames{strcmpi(colNames,'volume')} = 'volume';
        end
    end
    % Eliminate unused columns:
    for i = 1:length(data)
        colNames = data{i}.Properties.VariableNames;
        if sum(strcmpi(colNames,'AdjClose'))>0
            delInd = (strcmpi(colNames,'AdjClose') + strcmpi(colNames,'Date') + strcmpi(colNames,'volume'))==0;
            data{i}(:,delInd) = [];
            data{i} = data{i}(:,[find(strcmpi(data{i}.Properties.VariableNames,'Date')),...
                find(strcmpi(data{i}.Properties.VariableNames,'AdjClose')),...
                find(strcmpi(data{i}.Properties.VariableNames,'Volume'))]);
            data{i}.Properties.VariableNames{strcmpi(data{i}.Properties.VariableNames,'AdjClose')} = 'Close';
        elseif sum(strcmpi(colNames,'Close'))>0
            delInd = (strcmpi(colNames,'Close') + strcmpi(colNames,'Date') + strcmpi(colNames,'volume'))==0;
            data{i}(:,delInd) = [];
            data{i} = data{i}(:,[find(strcmpi(data{i}.Properties.VariableNames,'Date')),...
                find(strcmpi(data{i}.Properties.VariableNames,'Close')),...
                find(strcmpi(data{i}.Properties.VariableNames,'volume'))]);
        end
        if sum(strcmpi(data{i}.Properties.VariableNames,'volume'))
            data{i}.Properties.VariableNames{strcmpi(data{i}.Properties.VariableNames,'volume')} = 'volume';
        end
    end
    % Convert all dates to datenum's (i.e., number of days from January 0, 0000) (works irrespective of date format):
    sizes = zeros(1,length(data));
    for i = 1:length(data)
        data{i}.Date = datenum(data{i}.Date);
        sizes(i) = size(data{i},1);
    end
    % Find most limited data:
    latestStartingDate = 0;
    for i = 1:length(data)
        if min(data{i}.Date) > latestStartingDate
            latestStartingDate = min(data{i}.Date);
        end
    end
    % Eliminate data prior to this value for others:
    for i = 1:length(data)
        data{i}(data{i}.Date < latestStartingDate,:) = [];
    end
    % Find most limited data:
    earliestEndingDate = 9999999;
    for i = 1:length(data)
        if max(data{i}.Date) < earliestEndingDate
            earliestEndingDate = max(data{i}.Date);
        end
    end
    % Eliminate data prior to this value for others:
    for i = 1:length(data)
        data{i}(data{i}.Date > earliestEndingDate,:) = [];
    end
    % Check Rows for same dates:
    sizes = zeros(1,length(data));
    for i = 1:length(data)
        sizes(i) = size(data{i},1);
    end
    if ~all(sizes==mode(sizes))
       ixDifferent = find(sizes~=mode(sizes));
       warnStr = 'Some tickers (';
       for i = 1:length(ixDifferent)
           warnStr = [warnStr, tickers{i}, ', ']; %#ok
       end
       warnStr(end-1:end) = [];
       warnStr = [warnStr, ') do not have the same dates recorded. Trying interpolation.'];
       warning(warnStr);
        % If some dates different, try interpolating:
        ixUse = find(sizes==mode(sizes),1);
        for i = ixDifferent
            pp = csape(data{i}.Date, data{i}.Close, 'variational');
            pred = ppval(pp, data{ixUse}.Date); % evaluate PP at xx
            tempTable = data{ixUse}(:,1); tempTable.Close = pred;
            if isfield(data{i},'Volume')
                pp2 = csape(data{i}.Date, data{i}.Volume, 'variational');
                pred2 = ppval(pp2, data{ixUse}.Date); % evaluate PP at xx
                tempTable.volume = pred2;
            end
            data{i} = tempTable;
        end
    end
    % Convert cell2double, if any:
    for i = 1:length(data)
        if iscell(data{i}.Close)
            data{i}.Close = str2double(data{i}.Close);
        end
        if sum(strcmpi(data{i}.Properties.VariableNames,'volume')) && iscell(data{i}.volume)
            data{i}.volume = str2double(data{i}.volume);
        end
    end
    % Prepare output structure:
    out = struct('ticker',cell(1,length(data)),'date',cell(1,length(data)),...
        'raw',cell(1,length(data)),'return',cell(1,length(data)),'volume',...
        cell(1,length(data)));
    % Convert to return (daily percent change):
    for i = 1:length(data)
        returns = data{i}.Close(2:end,:) ./ data{i}.Close(1:end-1,:) - 1;
        out(i).ticker = tickers{i}; 
        out(i).date = data{i}.Date(2:end);
        out(i).raw = data{i}.Close(2:end); 
        out(i).return = returns;
        if sum(strcmpi(data{i}.Properties.VariableNames,'volume'))
            out(i).volume = data{i}.volume(2:end);
        end
    end
end