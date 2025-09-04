

% % ******************************************************
% % To get climatological and monthly 3D POC from cp660 profile data
% % Long-term Variability Analysis of Monthly POC Rates
% % ******************************************************

clear
clc

% Set file paths and parameters
filepath = '****\SCS_cp_profile\*****'; % 3D cp660 data path
Namelist = dir([filepath, '*.nc']);
LT = length(Namelist);
Sdep = [5, 10, 20, 30, 50, 75, 95, 115]; 
rx = 122;  % Grid dimensions for machine learning cp660 data:note this value should be changed when switch to optimazation derived cp660
ry = 199;  % Grid dimensions for machine learning cp660 data:note this value should be changed when switch to optimazation derived cp660

% Preallocate arrays for trend analysis
RCli = nan(rx, ry, length(Sdep));  % Trend slopes for each depth
PCli = nan(rx, ry, length(Sdep));  % p-values for each depth

% Preallocate arrays for climatological statistics
clmec_mean = nan(length(Sdep), LT);    % Mean POC for each depth and time
clmec_median = nan(length(Sdep), LT);  % Median POC for each depth and time
clmec_std = nan(length(Sdep), LT);     % Standard deviation for each depth and time

fprintf('Starting long-term variability analysis for %d depth levels...\n', length(Sdep));

% Process each depth level
for ds = 1:length(Sdep)
    fprintf('Processing depth level %d/%d: %.1f meters\n', ds, length(Sdep), Sdep(ds));
    
    % Initialize array for all time steps at current depth
    Alldata = nan(rx, ry, LT);
    
    %% Step 1: Extract data for current depth level
    for is = 1:LT
        % Read cp660 data
        da = ncread([filepath, Namelist(is).name], 'cp_profile_optcp4_NN');% here is the parameter for machine learing 3D cp660 files
        depth = ncread([filepath, Namelist(is).name], 'depth');
        
        % Find closest depth level
        [c, s] = min(abs(Sdep(ds) - depth));
        AS = da(:, :, s);
        AS(AS < 0) = nan;  % Remove negative values
        Alldata(:, :, is) = AS;
        
        % Display progress
        if mod(is, 50) == 0
            fprintf('  Depth %.1fm: Processed %d/%d files\n', Sdep(ds), is, LT);
        end
    end
    
    %% Step 2: Convert cp660 to POC concentration
    Alldata = 10.^((log10(Alldata) + 1.71)./1.27);
    Alldata1 = Alldata;  % Working copy for outlier removal
    
    %% Step 3: Outlier detection and removal (monthly basis)
    fprintf('  Performing outlier removal...\n');
    [ml, mj, ~] = size(Alldata1);
    
    for mot = 1:12
        % Create indices for all years of this month
        loc = mot:12:LT;
        
        AS1 = Alldata1(:, :, loc);
        
        for mmi = 1:ml
            for mmj = 1:mj
                BS = squeeze(AS1(mmi, mmj, :));
                
                if ~all(isnan(BS))
                    % Detect and remove outliers using median method
                    [BS1, TFrm, ~, ~, ~] = rmoutliers(BS, 'median');
                    
                    % Replace outliers with NaN if detected
                    if length(BS1) < length(BS)
                        BS(TFrm) = nan;
                        AS1(mmi, mmj, :) = BS;
                    end
                end
            end
        end
        
        % Update the original array
        Alldata1(:, :, loc) = AS1;
    end
    
    %% Step 4: Apply physical constraints and calculate statistics
    Alldata1(Alldata1 > 100/12) = nan;  % Remove unrealistically high values
    
    % Calculate temporal statistics for each time step
    for is = 1:LT
        btg = Alldata1(:, :, is) * 12;  % Convert to mg m-3 units 
        
        clmec_mean(ds, is) = mean(btg(:), 'omitnan');
        clmec_median(ds, is) = median(btg(:), 'omitnan');
        clmec_std(ds, is) = std(btg(:), 'omitnan');
    end
    
    %% Step 5: Calculate monthly climatology
    fprintf('  Calculating monthly climatology...\n');
    Cli_ave = nan(rx, ry, 12);
    
    for is = 1:12
        % Create indices for all years of this month
        loc1 = is:12:LT;
        Cli_ave(:, :, is) = mean(Alldata1(:, :, loc1), 3, 'omitnan');
    end
    
    %% Step 6: Calculate anomalies
    fprintf('  Calculating anomalies...\n');
    Ana_data = nan(rx, ry, LT);
    
    for dtg = 1:LT
        mtg = mod(dtg, 12);
        if (mtg == 0)
            mtg = 12;
        end
        Ana_data(:, :, dtg) = Alldata1(:, :, dtg) - Cli_ave(:, :, mtg);
    end
    
    %% Step 7: Calculate trends and significance
    fprintf('  Calculating trends and significance...\n');
    TNS = LT * 0.5;  % Minimum 50% data requirement
    Tml = 1:LT;      % Time vector
    
    for pd = 1:rx
        for pf = 1:ry
            % Extract time series for current grid cell
            t1 = Ana_data(pd, pf, :);
            t1 = squeeze(t1);
            
            % Remove NaN values
            valid_mask = ~isnan(t1);
            t12 = t1(valid_mask);
            t14 = Tml(valid_mask);
            
            % Check if sufficient data points are available
            if length(t12) > TNS
                % Prepare regression data
                X = [ones(length(t14), 1), t14'];
                Y = t12;
                
                % Perform linear regression
                [B, ~, ~, ~, ST] = regress(Y, X);
                
                % Store results
                RCli(pd, pf, ds) = B(2);  % Monthly trend
                PCli(pd, pf, ds) = ST(3); % p-value
            else
                % Insufficient data
                RCli(pd, pf, ds) = nan;
                PCli(pd, pf, ds) = nan;
            end
        end
        
        % Display progress
        if mod(pd, 20) == 0
            fprintf('    Row %d/%d completed\n', pd, rx);
        end
    end
    
    fprintf('Completed depth level %d/%d: %.1f meters\n', ds, length(Sdep), Sdep(ds));
end

%% Step 8: Additional Analysis and Output
fprintf('\nFinalizing analysis...\n');

% Calculate annual trends (monthly trend × 12)
RCli_annual = RCli * 12;

% Calculate significance masks (p < 0.05)
signif_mask = PCli < 0.05;

% Calculate percentage of significant trends for each depth
percent_signif = zeros(length(Sdep), 1);
for ds = 1:length(Sdep)
    valid_cells = sum(~isnan(PCli(:, :, ds)), 'all');
    if valid_cells > 0
        percent_signif(ds) = sum(signif_mask(:, :, ds), 'all') / valid_cells * 100;
    end
end

% Display summary statistics
fprintf('\n=== ANALYSIS SUMMARY ===\n');
fprintf('Depth levels analyzed: %d\n', length(Sdep));
fprintf('Time series length: %d months\n', LT);
fprintf('Minimum data requirement: %.0f%% of time series\n', 50);

fprintf('\nPercentage of significant trends (p < 0.05) by depth:\n');
for ds = 1:length(Sdep)
    fprintf('  %.1f m: %.1f%%\n', Sdep(ds), percent_signif(ds));
end

% Save results
save('SCS_POC_Monthly_Trend_Analysis.mat', ...
    'RCli', 'PCli', 'RCli_annual', 'signif_mask', ...
    'clmec_mean', 'clmec_median', 'clmec_std', ...
    'Sdep', 'LT', '-v7.3');

fprintf('\nAnalysis completed successfully!\n');
fprintf('Results saved to SCS_POC_Monthly_Trend_Analysis.mat\n');