clear
clc

% ******************************************************************
% Program to calculate climatological POC stock in the South China Sea
% 
% Three main phases:
% 1. Calculate monthly POC stock data
% 2. Handle outliers in monthly data
% 3. Compute climatological averages from corrected data
% ******************************************************************

%% Phase 1: Calculate Monthly POC Stock Data

% Depth levels for interpolation
isdepth = 5.01:0.5:130.01;

% Load bathymetry data for South China Sea study area
load('**\SCS\SCSdepth9km.mat');   
% Load predefined study region coordinates
load('**\SCS\SCSLatLon.mat'); 

% Define study region boundaries (Northern South China Sea)
[~, a1] = min(abs(lat - 15));
[~, b1] = min(abs(lat - 25.125));
[~, c1] = min(abs(lon - 105));
[~, d1] = min(abs(lon - 121.5));

% Set file paths
Cha_filepath = '****\Chla\****'; % Surface chlorophyll data path
filepath = '****\SCS_cp_profile\*****'; % 3D cp660 data path

% Preallocate arrays
Namelist = dir([filepath, '*', '.nc']);
LT = length(Namelist);
Clida = nan(122, 199, LT); % POC stock in euphotic zone
Clida_mean = nan(122, 199, LT); % Mean POC in euphotic zone
Clida_100 = nan(122, 199, LT); % POC stock in 100m layer
Clida_100_mean = nan(122, 199, LT); % Mean POC in 100m layer

% Process each time step
for is = 1:LT
    yg = Namelist(is).name;
    
    % Read 3D cp660 data (machine learning variable)
    da = ncread([filepath, yg], 'cp_profile_optcp4_NN');  % here is the parameter for machine learing 3D cp660 files
    depth = ncread([filepath, yg], 'depth');
    
    % Read surface chlorophyll data
    ChlaNaList = dir([Cha_filepath, '*MODIS.', yg(19:24), '*9km.nc']);  % get the month information 
    if ~isempty(ChlaNaList)
        ChlaName = ChlaNaList(1).name;
        Chla = ncread([Cha_filepath, ChlaName], 'chlor_a', [3395, 779], [232, 288], [1, 1])'; % get data for the whole SCS
        Chla = Chla(b1:a1, c1:d1); % Extract Northern SCS region
        
        % Calculate euphotic zone depth using Shang et al.(2011) model
        Zeu2 = 10.^(1.35 - 0.4026.*log10(Chla) + 0.0375.*log10(Chla).^2);
        
        % Handle negative values
        Zeu2(Zeu2 < 0) = nan;
        da(da < 0) = nan;
        
        % Process each grid cell
        [ml, mj, ~] = size(da);
        for mmi = 1:ml
            for mmj = 1:mj
                if ~isnan(Chla(mmi, mmj)) && ~isnan(Zeu2(mmi, mmj))
                    % Extract cp660 profile
                    BS = da(mmi, mmj, :);
                    BS = squeeze(BS);
                    
                    % Convert cp660 to POC concentration (umol/L)
                    isPOC = 10.^((log10(BS) + 1.71)./1.27);
                    
                    % Define integration depths
                    StZeu = 1:1:min(Zeu2(mmi, mmj), max(depth));
                    StZeu2 = 1:1:min(100, max(depth));
                    
                    % Interpolate and integrate POC profiles
                    if length(depth) > 1 && ~all(isnan(isPOC))
                        % Euphotic zone integration
                        interpPOC = interp1([1; depth], [isPOC(1); isPOC], StZeu, 'linear', 'extrap');
                        Clida(mmi, mmj, is) = sum(interpPOC, 'omitnan');
                        Clida_mean(mmi, mmj, is) = mean(interpPOC, 'omitnan');
                        
                        % 100m layer integration
                        interpPOC_100 = interp1([1; depth], [isPOC(1); isPOC], StZeu2, 'linear', 'extrap');
                        Clida_100(mmi, mmj, is) = sum(interpPOC_100, 'omitnan');
                        Clida_100_mean(mmi, mmj, is) = mean(interpPOC_100, 'omitnan');
                    end
                end
            end
        end
    end
    
    % Display progress
    fprintf('Processed file %d of %d: %s\n', is, LT, yg);
end

%% Phase 2: Outlier Detection and Removal in Monthly POC Data

% Process each month separately
for mot = 1:12
    % Create indices for all years of this month
    loc = mot:12:(12*21);
    
    % Extract data for this month across all years
    AS1 = Clida(:, :, loc); % Euphotic zone stock
    AS3 = Clida_mean(:, :, loc); % Euphotic zone mean
    AS2 = Clida_100(:, :, loc); % 100m layer stock
    AS5 = Clida_100_mean(:, :, loc); % 100m layer mean
    
    % Apply physical constraints
    AS1(AS1 < 0) = nan;
    AS3(AS3 < 0) = nan;
    AS1(AS1 > 5000/12) = nan; % Reasonable upper limit
    AS3(AS3 > 100/12) = nan;
    
    AS2(AS2 < 0) = nan;
    AS5(AS5 < 0) = nan;
    AS2(AS2 > 5000/12) = nan;
    AS5(AS5 > 100/12) = nan;
    
    % Process each grid cell for outlier removal
    [ml, mj, ~] = size(AS1);
    for mmi = 1:ml
        for mmj = 1:mj
            % Process euphotic zone stock
            BS = squeeze(AS1(mmi, mmj, :));
            if ~all(isnan(BS))
                [BS1, TFrm] = rmoutliers(BS, 'median');
                if length(BS1) < length(BS) % If outliers were detected
                    BS(TFrm) = nan;
                end
                AS1(mmi, mmj, :) = BS;
            end
            
            % Process euphotic zone mean
            BS = squeeze(AS3(mmi, mmj, :));
            if ~all(isnan(BS))
                [BS1, TFrm] = rmoutliers(BS, 'median');
                if length(BS1) < length(BS)
                    BS(TFrm) = nan;
                end
                AS3(mmi, mmj, :) = BS;
            end
            
            % Process 100m layer stock
            BS = squeeze(AS2(mmi, mmj, :));
            if ~all(isnan(BS))
                [BS1, TFrm] = rmoutliers(BS, 'median');
                if length(BS1) < length(BS)
                    BS(TFrm) = nan;
                end
                AS2(mmi, mmj, :) = BS;
            end
            
            % Process 100m layer mean
            BS = squeeze(AS5(mmi, mmj, :));
            if ~all(isnan(BS))
                [BS1, TFrm] = rmoutliers(BS, 'median');
                if length(BS1) < length(BS)
                    BS(TFrm) = nan;
                end
                AS5(mmi, mmj, :) = BS;
            end
        end
    end
    
    % Store processed data back
    Clida(:, :, loc) = AS1;
    Clida_mean(:, :, loc) = AS3;
    Clida_100(:, :, loc) = AS2;
    Clida_100_mean(:, :, loc) = AS5;
    
    fprintf('Completed outlier removal for month %d\n', mot);
end

%% Phase 3: Compute Climatological Averages

% Preallocate climatology arrays
climato1 = nan(ml, mj, 12); % Euphotic zone stock climatology
climato2 = nan(ml, mj, 12); % 100m layer stock climatology
climato3 = nan(ml, mj, 12); % Euphotic zone mean climatology
climato4 = nan(ml, mj, 12); % 100m layer mean climatology

% Calculate monthly climatologies
for mot = 1:12
    % Create indices for all years of this month
    loc = mot:12:(12*21);
    
    % Compute climatological means
    climato1(:, :, mot) = mean(Clida(:, :, loc), 3, 'omitnan');
    climato2(:, :, mot) = mean(Clida_100(:, :, loc), 3, 'omitnan');
    climato3(:, :, mot) = mean(Clida_mean(:, :, loc), 3, 'omitnan');
    climato4(:, :, mot) = mean(Clida_100_mean(:, :, loc), 3, 'omitnan');
    
    fprintf('Computed climatology for month %d\n', mot);
end

disp('Climatological POC calculation completed successfully');


%% Phase 4: Calculate Long-term Trend and Statistical Significance

% Get dimensions of climatology data
[rx, ry, ~] = size(climato1);

% Calculate anomalies by removing climatological seasonal cycle
Ana_data = nan(rx, ry, 252);       % Anomalies for euphotic zone stock
Ana_100_data = nan(rx, ry, 252);   % Anomalies for 100m layer stock

fprintf('Calculating anomalies...\n');
for dtg = 1:252
    % Determine corresponding month (1-12)
    mtg = mod(dtg, 12);
    if (mtg == 0)
        mtg = 12;
    end
    
    % Calculate anomalies: observed minus climatology
    Ana_data(:, :, dtg) = Clida(:, :, dtg) - climato1(:, :, mtg);
    Ana_100_data(:, :, dtg) = Clida_100(:, :, dtg) - climato2(:, :, mtg);
    
    % Display progress every 50 time steps
    if mod(dtg, 50) == 0
        fprintf('Anomaly calculation: %d/%d completed\n', dtg, 252);
    end
end

% Initialize arrays for trend analysis
RCli = nan(rx, ry);        % Trend slope for euphotic zone
RCli_100 = nan(rx, ry);    % Trend slope for 100m layer
PCli = nan(rx, ry);        % p-value for euphotic zone trend
PCli_100 = nan(rx, ry);    % p-value for 100m layer trend

% Set minimum data requirement threshold (50% of total time series)
TNS = 252 * 0.5;  % Minimum 126 valid time points required
Tml = 1:252;      % Time vector (months)

fprintf('Calculating trends and significance levels...\n');

% Process each grid cell
for pd = 1:rx
    for pf = 1:ry
        %% Euphotic Zone Trend Analysis
        % Extract time series for current grid cell
        t1 = Ana_data(pd, pf, :);
        t1 = squeeze(t1);
        
        % Remove NaN values
        valid_mask = ~isnan(t1);
        t12 = t1(valid_mask);
        t14 = Tml(valid_mask);  % Corresponding time points
        
        % Check if sufficient data points are available
        if length(t12) > TNS
            % Prepare regression data
            X = [ones(length(t14), 1), t14'];  % Design matrix [intercept, time]
            Y = t12;                           % Response variable (anomalies)
            
            % Perform linear regression
            [B, ~, ~, ~, ST] = regress(Y, X);
            
            % Store results
            RCli(pd, pf) = B(2);  % Slope (trend per month)
            PCli(pd, pf) = ST(3); % p-value for slope
            
            % Alternative: Use polyfit for additional statistics
            % [p, S] = polyfit(t14, t12, 1);
            % RCli(pd, pf) = p(1);
            % [~, PCli(pd, pf)] = corr(t14', t12);
        else
            % Insufficient data
            RCli(pd, pf) = nan;
            PCli(pd, pf) = nan;
        end
        
        %% 100m Layer Trend Analysis
        % Extract time series for current grid cell
        t1_100 = Ana_100_data(pd, pf, :);
        t1_100 = squeeze(t1_100);
        
        % Remove NaN values
        valid_mask_100 = ~isnan(t1_100);
        t12_100 = t1_100(valid_mask_100);
        t14_100 = Tml(valid_mask_100);
        
        % Check if sufficient data points are available
        if length(t12_100) > TNS
            % Prepare regression data
            X_100 = [ones(length(t14_100), 1), t14_100'];
            Y_100 = t12_100;
            
            % Perform linear regression
            [B_100, ~, ~, ~, ST_100] = regress(Y_100, X_100);
            
            % Store results
            RCli_100(pd, pf) = B_100(2);    % Slope
            PCli_100(pd, pf) = ST_100(3);   % p-value
        else
            % Insufficient data
            RCli_100(pd, pf) = nan;
            PCli_100(pd, pf) = nan;
        end
    end
    
    % Display progress every 10 rows
    if mod(pd, 10) == 0
        fprintf('Trend analysis: row %d/%d completed\n', pd, rx);
    end
end

%% Additional Statistical Analysis and Output

% Convert monthly trends to annual trends (multiply by 12 months)
RCli_annual = RCli * 12;
RCli_100_annual = RCli_100 * 12;

% Calculate significance masks (p < 0.05)
signif_mask = PCli < 0.05;
signif_mask_100 = PCli_100 < 0.05;

% Calculate percentage of significant trends
percent_signif = sum(signif_mask(:), 'omitnan') / sum(~isnan(PCli(:))) * 100;
percent_signif_100 = sum(signif_mask_100(:), 'omitnan') / sum(~isnan(PCli_100(:))) * 100;

fprintf('Trend analysis completed:\n');
fprintf('Euphotic zone: %.1f%% of grid cells show significant trends (p < 0.05)\n', percent_signif);
fprintf('100m layer: %.1f%% of grid cells show significant trends (p < 0.05)\n', percent_signif_100);

% Calculate basic statistics for trends
fprintf('\nEuphotic Zone Trend Statistics:\n');
fprintf('Mean trend: %.4f units/month (%.4f units/year)\n', ...
    mean(RCli(:), 'omitnan'), mean(RCli_annual(:), 'omitnan'));
fprintf('Trend range: [%.4f, %.4f] units/month\n', ...
    min(RCli(:), 'omitnan'), max(RCli(:), 'omitnan'));

fprintf('\n100m Layer Trend Statistics:\n');
fprintf('Mean trend: %.4f units/month (%.4f units/year)\n', ...
    mean(RCli_100(:), 'omitnan'), mean(RCli_100_annual(:), 'omitnan'));
fprintf('Trend range: [%.4f, %.4f] units/month\n', ...
    min(RCli_100(:), 'omitnan'), max(RCli_100(:), 'omitnan'));

%% Optional: Save Results
% save('SCS_POC_Trend_Analysis.mat', ...
%     'RCli', 'PCli', 'RCli_100', 'PCli_100', ...
%     'RCli_annual', 'RCli_100_annual', ...
%     'signif_mask', 'signif_mask_100', ...
%     'Ana_data', 'Ana_100_data', ...
%     'lat', 'lon');

disp('Phase 4: Long-term trend analysis completed successfully');

