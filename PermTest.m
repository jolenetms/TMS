clear all
close all

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 31, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["StudyID", "CONDITION", "HRSD_Tot_W0", "HRSD_Tot_W1", "HRSD_Tot_W2", "HRSD_Tot_W3", "HRSD_Tot_W4", "HRSD_Tot_W5", "HRSD_Tot_Final", "HRSD_Tot_FU1", "HRSD_Tot_FU2", "HRSD_Tot_FU3", "IDS_Tot_W0", "IDS_Tot_W1", "IDS_Tot_W2", "IDS_Tot_W3", "IDS_Tot_W4", "IDS_Tot_W5", "IDS_Tot_Final", "IDS_Tot_FU1", "IDS_Tot_FU2", "IDS_Tot_FU3", "QIDS_Tot_W0", "QIDS_Tot_W1", "QIDS_Tot_W2", "QIDS_Tot_W3", "QIDS_Tot_W4", "QIDS_Tot_W5", "QIDS_Tot_Final", "QIDS_Tot_FU1", "QIDS_Tot_FU2", "QIDS_Tot_FU3", "BENZO"];
opts.VariableTypes = ["string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "StudyID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "StudyID", "EmptyFieldRule", "auto");

% Import the data
TESTS_by_week = readtable("/Users/jolinchou/Documents/Sudden_Gains/Three-D/IndiraData0129.csv", opts);
TESTS_by_week_array = table2array(TESTS_by_week(:, 3:end));

%% Clear temporary variables
clear opts

%% Initalize Variables

% Initialize a logical array to store results
meetsCriteria = false(size(TESTS_by_week_array(:, 2:8))); % Assuming 372 rows
 allHAMD_Delta = zeros(size(TESTS_by_week_array)); 
 sg_reversalvalue = [];

% Initialize logical array to store remission data starting at week 2
remission_criteria_met = false(size(TESTS_by_week_array, 1), 9);
response_criteria_met = false(size(TESTS_by_week_array, 1), 10);

%% Detect remission 
for p = 1:size(TESTS_by_week_array, 1)
    participant_data = TESTS_by_week_array(p, 1:10);
    
    % Create value to check if two consecutive HAMD values are <= 10, and that any consecutive increase is either <= 3, or value of <= 6
    remission_criteria_met_10 = (participant_data(1:end-1) <= 10) & (participant_data(2:end) <= 10);
    increase = participant_data(2:end) - participant_data(1:end-1);
    second_less_than_6 = participant_data(2:end) <= 6;
    remission_criteria_met_ppt = remission_criteria_met_10 & ((increase <= 3) | second_less_than_6);
    
    % Create variables for a HAMD of half baseline HAMD to compare to the values at the following weeks
    halfbaseline = participant_data(:, 1) ./ 2;
    response_criteria = participant_data <= halfbaseline;
    
    % Check all remission criteria for each participant
    remission_criteria_met(p, :) = remission_criteria_met_ppt;
    response_criteria_met(p, :) = response_criteria;
end

%% Normalize data table's number of weeks by adding columns between follow ups 

% Normalize FU1 to FU2
[m, n] = size(remission_criteria_met);

% Specify the index where you want to add three columns (between columns 7 and 8)
insertIndex = 7;

% Initialize a new matrix with three additional columns
newMatrix = zeros(m, n + 2);

% Iterate through each row
for row = 1:m
    % Copy the columns before the insertion point
    newMatrix(row, 1:insertIndex) = remission_criteria_met(row, 1:insertIndex);
    newMatrix(row, (insertIndex + 3):end) = remission_criteria_met(row, (insertIndex + 1):end);
    
    % Check the condition for the surrounding columns
    if insertIndex > 1 && insertIndex < n && remission_criteria_met(row, insertIndex) == 1 && remission_criteria_met(row, insertIndex + 1) == 1
        % If 1's on both sides, add three columns of 1's
        newMatrix(row, insertIndex+1:insertIndex+2) = 1;
    else
        % Otherwise, add three columns of 0's
        newMatrix(row, insertIndex+1:insertIndex+2) = 0;
    end
    
    % Copy the columns after the insertion point
    newMatrix(row, insertIndex+3:end) = remission_criteria_met(row, insertIndex+1:end);
end

%Normalize FU2 to FU3
[o, p] = size(newMatrix);

% Specify the index where you want to add three columns (between columns 7 and 8)
insertIndex2 = 10;

% Initialize a new matrix with three additional columns
final_remission_crit_met = zeros(o, p + 7);

% Iterate through each row
for row = 1:o
    % Copy the columns before the insertion point
    final_remission_crit_met(row, 1:insertIndex2) = newMatrix(row, 1:insertIndex2);
    final_remission_crit_met(row, insertIndex2 + 8:end) = newMatrix(row, insertIndex2 + 1:end);
    
    % Check the condition for the surrounding columns
    if insertIndex2 > 1 && insertIndex2 < p && newMatrix(row, insertIndex2) == 1 && newMatrix(row, insertIndex2 + 1) == 1
        % If 1's on both sides, add three columns of 1's
        final_remission_crit_met(row, (insertIndex2 + 1): (insertIndex2 + 7)) = 1;
    else
        % Otherwise, add three columns of 0's
        final_remission_crit_met(row, (insertIndex2 + 1) : (insertIndex2 + 7)) = 0;
    end
    
    % Copy the columns after the insertion point
    final_remission_crit_met(row, (insertIndex2 + 8):end) = newMatrix(row, (insertIndex2 + 1):end);
end

%% Find remitters with sustained remission

twowk_sustainedremission = sum(conv2(final_remission_crit_met, [1 1 1], 'valid') == 3, 2) > 0;
twowk_ppt_rem = find(twowk_sustainedremission);

fourwk_sustainedremission = sum(conv2(final_remission_crit_met, [1 1 1 1 1], 'valid') == 5, 2) > 0;
fourwk_ppt_rem = find(fourwk_sustainedremission);

sixwk_sustainedremission = sum(conv2(final_remission_crit_met, [1 1 1 1 1 1 1], 'valid') == 7, 2) > 0;
sixwk_ppt_rem = find(sixwk_sustainedremission);

eightwk_sustainedremission = sum(conv2(final_remission_crit_met, [1 1 1 1 1 1 1 1 1], 'valid') == 9, 2) > 0;
eightwk_ppt_rem = find(eightwk_sustainedremission);

twelvewk_sustainedremission = sum(conv2(final_remission_crit_met, [1 1 1 1 1 1 1 1 1 1 1 1 1], 'valid') == 13, 2) > 0;
twelvewk_ppt_rem = find(twelvewk_sustainedremission);

%% General Analysis of Remission & Response

response_criteria_met_atall = false(size(remission_criteria_met, 1), 1); % initialize variable

% Check if remission & response criteria is met at least once for each ppt
remission_criteria_met_atall = any(remission_criteria_met == 1, 2);
response_criteria_met_atall = any(response_criteria_met == 1, 2);

% Get sum of total remitters & responders
total_remitters = sum(remission_criteria_met_atall);
total_responders = sum(response_criteria_met_atall);

%% Index each participants first instance of remission

% Extract the row (ppt) and column (week) that remission & response was reached for each instance of remission
[row_remission_index, col_remission_index] = find(remission_criteria_met == 1);
[row_response_index, col_response_index] = find(response_criteria_met == 1);

% Correct the session data to align with the original data by adding 1 the indexed sessions
% col_remission_index_corrected = col_remission_index - 1;

% Concatenize the ppt and session data into the same 2-column variable and sort it so that it is in order through the ppts
remission_instances = horzcat(row_remission_index, col_remission_index);
remission_instances_sorted = sortrows(remission_instances, 1);
response_instances = horzcat(row_response_index, col_response_index);
response_instances_sorted = sortrows(response_instances);

% Extract the first instances of remission & response for each ppt
[first_rem_instance, rem_weekidx] = unique(remission_instances_sorted(:, 1), 'stable');
first_rem_week = remission_instances_sorted(rem_weekidx, 2);
[first_resp_instance, resp_weekidx] = unique(response_instances_sorted(:, 1), 'stable');
first_resp_week = response_instances_sorted(resp_weekidx, 2);
% Concatenize first remission instance data into one 2-column variable
remission_met_indices = horzcat(first_rem_instance, first_rem_week);
response_met_indices = horzcat(first_resp_instance, first_resp_week);

remission_criteria_week = zeros(size(remission_criteria_met, 1), 1);
response_criteria_week = zeros(size(response_criteria_met, 1), 1);

remission_criteria_week(first_rem_instance) = first_rem_week;
response_criteria_week(first_resp_instance) = first_resp_week;

%% Evaluate sudden gains for each ppt
% Loop through each participant
for weekindex = 2:5
    % Compute differences between consecutive scores for each participant
    HAMD_Delta =  TESTS_by_week_array(:, weekindex + 1) - TESTS_by_week_array(:, weekindex);

    allHAMD_Delta(:, weekindex) = HAMD_Delta;

    halfdelta = allHAMD_Delta./2; %divide the change in HAMD by 2

    reversalvalue = TESTS_by_week_array + halfdelta;

    % Calculate the percentage differences
    HAMD_percentdecrease =  (TESTS_by_week_array(:, weekindex) - TESTS_by_week_array(:, weekindex + 1)) ./ TESTS_by_week_array(:, weekindex) * 100;

    % Extract 3 Pre and 3 Post HAMD values from gain
    if weekindex == 2
         HAMD_pre = TESTS_by_week_array(:, (weekindex - 1):(weekindex)); 
    else 
        HAMD_pre = TESTS_by_week_array(:, (weekindex - 2):(weekindex)); 
    end

    if weekindex == 5
         HAMD_post = TESTS_by_week_array(:, (weekindex + 1):(weekindex + 2)); 
    else 
        HAMD_post = TESTS_by_week_array(:, (weekindex + 1):(weekindex + 3));
    end  

    % Calculate pre standard deviation 
    std_pre = std((HAMD_pre), 0, 2, 'omitmissing');

    % Calculate post standard deviation 
    std_post = std((HAMD_post), 0, 2, 'omitmissing');

    % Assuming HAMD_pre and HAMD_post are your arrays
    rows_with_NaN = any(isnan(HAMD_pre), 2) | any(isnan(HAMD_post), 2);

    % Assuming HAMD_pre and HAMD_post are your arrays
    rows_with_NaN_pre = sum(isnan(HAMD_pre), 2);
    rows_with_NaN_post = sum(isnan(HAMD_post), 2);

    % Rows with 0 NaN
    rows_with_zeros = (rows_with_NaN_pre == 0) & (rows_with_NaN_post == 0);

    % Rows with 1 NaN
    rows_with_1_NaN = (rows_with_NaN_pre == 1) & (rows_with_NaN_post == 0) | (rows_with_NaN_pre == 0) & (rows_with_NaN_post == 1);

    % Rows with 1 NaN pre and post
    rows_with_1_NaN_prepost = (rows_with_NaN_pre == 1) & (rows_with_NaN_post == 1);

    % Rows with 2 NaNs pre or post
    rows_with_2_NaNs = (rows_with_NaN_pre >= 2) | (rows_with_NaN_post >= 2);

    % Indexing with rows_with_NaN to access rows with NaN values
    zero_rows = find(rows_with_zeros);
    One_NaN_rows = find(rows_with_1_NaN);
    PrePost_NaN_rows = find(rows_with_1_NaN_prepost);
    bye_rows =  find(rows_with_2_NaNs);

     %% Double Check
  
    % Create a range of all possible row indexes
    totalRows = 1:372; % Replace 'totalNumberOfRows' with the total number of rows in your data

     % Combine all indexes into one array
    allIndexes = sort([zero_rows; One_NaN_rows; PrePost_NaN_rows; bye_rows]);

    % Find the missing rows
    missingRows = setdiff(totalRows, unique(allIndexes));

    doublecheck = length(One_NaN_rows) + length(PrePost_NaN_rows) + length(bye_rows) + length(zero_rows);

    %% Apply criteria 3, adjusting for missing values
            
    allsymptomfluctuation(bye_rows,:) = NaN;

    % Apply symptom fluctuation formula for one missing value pre and post
        for i = 1:length(PrePost_NaN_rows)
            PrePost_index = PrePost_NaN_rows(i);
            symptomfluctuation = 4.303 * sqrt(((TESTS_by_week_array(PrePost_index, weekindex) - 1) .* std_pre(PrePost_index).^2 + (TESTS_by_week_array(PrePost_index, weekindex + 1) - 1) .* std_post(PrePost_index).^2) ./ (TESTS_by_week_array(PrePost_index, weekindex) + TESTS_by_week_array(PrePost_index, weekindex + 1) - 2));

            allsymptomfluctuation(PrePost_index,:) = symptomfluctuation;
        end 

        % Apply symptom fluctuation formula for one missing value 
        for j = 1:length(One_NaN_rows)
            OneNaN_index = One_NaN_rows(j);
            symptomfluctuation = 3.182 * sqrt(((TESTS_by_week_array(OneNaN_index, weekindex) - 1) .* std_pre(OneNaN_index).^2 + (TESTS_by_week_array(OneNaN_index, weekindex + 1) - 1) .* std_post(OneNaN_index).^2) ./ (TESTS_by_week_array(OneNaN_index, weekindex) + TESTS_by_week_array(OneNaN_index, weekindex + 1) - 2));

            allsymptomfluctuation(OneNaN_index,:) = symptomfluctuation;
        end 

        % Apply symptom fluctuation formula for no missing values
        for k = 1:length(zero_rows)
            zero_index = zero_rows(k);
        symptomfluctuation = 2.776 * sqrt(((TESTS_by_week_array(zero_index, weekindex) - 1) .* std_pre(zero_index).^2 + (TESTS_by_week_array(zero_index, weekindex + 1) - 1) .* std_post(zero_index).^2) ./ (TESTS_by_week_array(zero_index, weekindex) + TESTS_by_week_array(zero_index, weekindex + 1) - 2));

         allsymptomfluctuation(zero_index,:) = symptomfluctuation;
        end

         % Criteria 3: Calculate Mpre - Mpost
    MeanDelta = mean(HAMD_pre, 2, 'omitmissing') - mean(HAMD_post, 2, 'omitmissing');
    
    %Check if all three criteria are met for each score difference
    criteriaCheck = HAMD_Delta <= -7 & HAMD_percentdecrease >= 25 &  MeanDelta > allsymptomfluctuation;

    meetsCriteria(:, weekindex) = criteriaCheck;
  
end

%% Find reversals occuring for each participant

% Index the sudden gain events
[rowIndex, colIndex] = find(meetsCriteria == 1);
sg_linearindex = sub2ind(size(meetsCriteria), rowIndex, colIndex);

% Find reversal value corresponding to the sudden gain index
sg_reversalvalue = reversalvalue(sg_linearindex);

% Concatenate sudden gain index with reversal value
sg_index = horzcat(rowIndex, colIndex, sg_reversalvalue);

% Find indices of sessions that occur after the sudden gain/find indices of zeros occurring after ones for each row
zero_after_one_indices = cell(size(meetsCriteria, 1), 1);
for row = 1:size(meetsCriteria, 1)
    % Find where ones occur
    one_indices = find(meetsCriteria(row, :) == 1);

    % Initialize array to store indices of zeros after ones
    zero_indices = [];

    % Collect indices of zeros after each one
    for i = 1:length(one_indices)
        if one_indices(i) < length(meetsCriteria(row, :))
            zero_indices = [zero_indices, find(meetsCriteria(row, one_indices(i):end) == 0) + one_indices(i) - 1];
        end
    end

    zero_after_one_indices{row} = zero_indices; 

% Remove duplicates and store in cell array
    zero_after_one_indices{row} = unique(zero_indices); % Find index of values after the sudden gain
end

%initalize variable to store matrix with all reversals
consolidatedreversal = false(size(TESTS_by_week_array));

% Compare HAMD values collected after the sudden gain to the reversal value
% rowIndex indicates rows where sudden gain occured
for i = 1:numel(rowIndex)

    % Store indices of values after the sudden gain
    indices = zero_after_one_indices{rowIndex(i)};
    
    % sg_index(i,2) = column index where sudden gain even toccured
    indices = indices(indices >= sg_index(i,2));

    % If the HAMD is greater than reversal value, store in reversal_ppt
reversal_ppt = TESTS_by_week_array(rowIndex(i), indices) >= sg_index(i,3);

consolidatedreversal(rowIndex(i), indices) = reversal_ppt;
end 

rows_with_ones_sg = any(meetsCriteria == 1, 2);
rows_with_ones_reversals = any(consolidatedreversal == 1, 2);

indices_with_ones = find(rows_with_ones_sg);
indices_with_ones_reversals = find(rows_with_ones_reversals);

% Create a new variable and assign modified values
rows_with_ones_sg_reversed = rows_with_ones_sg;
rows_with_ones_sg_reversed(indices_with_ones_reversals) = 0; % Sudden gainers accounting for reversal criteria

%% Permutation Test

% Create remitter variable to sample from for permuation test
remitters_sg_nonsg = rows_with_ones_sg_reversed(remission_criteria_met_atall, :); % total sample

twowk_sus_rem_sg = sum(rows_with_ones_sg_reversed(twowk_sustainedremission));
twowk_sus_rem_sg_nonsg = sum(twowk_sustainedremission);

twowk_sus_rem_sg_nonsg = rows_with_ones_sg_reversed(twowk_sustainedremission, :);
twowk_sus_rem_observedratio = mean(twowk_sus_rem_sg_nonsg);

fourwk_sus_rem_sg_nonsg = rows_with_ones_sg_reversed(fourwk_sustainedremission, :);
fourwk_sus_rem_observedratio = mean(fourwk_sus_rem_sg_nonsg);

sixwk_sus_rem_sg_nonsg = rows_with_ones_sg_reversed(sixwk_sustainedremission, :);
sixwk_sus_rem_observedratio = mean(sixwk_sus_rem_sg_nonsg);

eightwk_sus_rem_sg_nonsg = rows_with_ones_sg_reversed(eightwk_sustainedremission, :);
eightwk_sus_rem_observedratio = mean(eightwk_sus_rem_sg_nonsg);

twelvewk_sus_rem_sg_nonsg = rows_with_ones_sg_reversed(twelvewk_sustainedremission, :);
twelvewk_sus_rem_observedratio = mean(twelvewk_sus_rem_sg_nonsg);

% Number of permutations
num_permutations = 1000000;

% Initialize array to store permuted test statistics
permuted_statistics = zeros(num_permutations, 1);

% Perform permutation testing
for i = 1:num_permutations
    % Permute the data
    shuffled_data = remitters_sg_nonsg(randperm(154, 100));
    %x(:,i) = shuffled_data

    % Compute test statistic for permuted data (mean difference)
    permuted_statistics(i) = mean(shuffled_data);
end

% Calculate p-value
p_value_twowk = sum(permuted_statistics >= twowk_sus_rem_observedratio) / num_permutations;
   %p_value2 = mean(permuted_statistics >= twowk_sus_rem_observedratio);
p_value_fourwk = sum(permuted_statistics >= fourwk_sus_rem_observedratio) / num_permutations;
p_value_sixwk = sum(permuted_statistics >= sixwk_sus_rem_observedratio) / num_permutations;
p_value_eightwk = sum(permuted_statistics >= eightwk_sus_rem_observedratio) / num_permutations;
p_value_twelvewk = sum(permuted_statistics >= twelvewk_sus_rem_observedratio) / num_permutations;

%% 2 week sustained remission p value
if twowk_sus_rem_observedratio < median(permuted_statistics)
    direction = 'left';  % Left tail test statistic
elseif twowk_sus_rem_observedratio > median(permuted_statistics)
    direction = 'right';  % Right tail test statistic
else
    direction = 'both';  % Test statistic at center (e.g., mean difference is exactly 0)
end

% Compute p-value based on the direction of the observed test statistic
if strcmp(direction, 'left')
    p_value_twowk = sum(permuted_statistics <= twowk_sus_rem_observedratio) / num_permutations;
elseif strcmp(direction, 'right')
    p_value_twowk = sum(permuted_statistics >= twowk_sus_rem_observedratio) / num_permutations;
else
    p_value_twowk = 2 * min(sum(permuted_statistics <= twowk_sus_rem_observedratio), sum(permuted_statistics >= twowk_sus_rem_observedratio)) / num_permutations;
end

%% 4 week sustained remission p value
% Determine the direction of the observed test statistic
if fourwk_sus_rem_observedratio < median(permuted_statistics)
    direction = 'left';  % Left tail test statistic
elseif fourwk_sus_rem_observedratio > median(permuted_statistics)
    direction = 'right';  % Right tail test statistic
else
    direction = 'both';  % Test statistic at center (e.g., mean difference is exactly 0)
end

% Compute p-value based on the direction of the observed test statistic
if strcmp(direction, 'left')
    p_value_fourwk = sum(permuted_statistics <= fourwk_sus_rem_observedratio) / num_permutations;
elseif strcmp(direction, 'right')
    p_value_fourwk = sum(permuted_statistics >= fourwk_sus_rem_observedratio) / num_permutations;
else
    p_value_fourwk = 2 * min(sum(permuted_statistics <= fourwk_sus_rem_observedratio), sum(permuted_statistics >= fourwk_sus_rem_observedratio)) / num_permutations;
end


%% 6 week sustained remission p value
% Determine the direction of the observed test statistic
if sixwk_sus_rem_observedratio < median(permuted_statistics)
    direction = 'left';  % Left tail test statistic
elseif sixwk_sus_rem_observedratio > median(permuted_statistics)
    direction = 'right';  % Right tail test statistic
else
    direction = 'both';  % Test statistic at center (e.g., mean difference is exactly 0)
end

% Compute p-value based on the direction of the observed test statistic
if strcmp(direction, 'left')
    p_value_sixwk = sum(permuted_statistics <= sixwk_sus_rem_observedratio) / num_permutations;
elseif strcmp(direction, 'right')
    p_value_sixwk = sum(permuted_statistics >= sixwk_sus_rem_observedratio) / num_permutations;
else
    p_value_sixwk = 2 * min(sum(permuted_statistics <= sixwk_sus_rem_observedratio), sum(permuted_statistics >= sixwk_sus_rem_observedratio)) / num_permutations;
end

%% 8 week sustained remission p value
% Determine the direction of the observed test statistic
if eightwk_sus_rem_observedratio < median(permuted_statistics)
    direction = 'left';  % Left tail test statistic
elseif eightwk_sus_rem_observedratio > median(permuted_statistics)
    direction = 'right';  % Right tail test statistic
else
    direction = 'both';  % Test statistic at center (e.g., mean difference is exactly 0)
end

% Compute p-value based on the direction of the observed test statistic
if strcmp(direction, 'left')
    p_value_eightwk = sum(permuted_statistics <= eightwk_sus_rem_observedratio) / num_permutations;
elseif strcmp(direction, 'right')
    p_value_eightwk = sum(permuted_statistics >= eightwk_sus_rem_observedratio) / num_permutations;
else
    p_value_eightwk = 2 * min(sum(permuted_statistics <= eightwk_sus_rem_observedratio), sum(permuted_statistics >= eightwk_sus_rem_observedratio)) / num_permutations;
end

%% 12 week sustained remission p value
% Determine the direction of the observed test statistic
if twelvewk_sus_rem_observedratio < median(permuted_statistics)
    direction = 'left';  % Left tail test statistic
elseif twelvewk_sus_rem_observedratio > median(permuted_statistics)
    direction = 'right';  % Right tail test statistic
else
    direction = 'both';  % Test statistic at center (e.g., mean difference is exactly 0)
end

% Compute p-value based on the direction of the observed test statistic
if strcmp(direction, 'left')
    p_value_twelvewk = sum(permuted_statistics <= twelvewk_sus_rem_observedratio) / num_permutations;
elseif strcmp(direction, 'right')
    p_value_twelvewk = sum(permuted_statistics >= twelvewk_sus_rem_observedratio) / num_permutations;
else
    p_value_twelvewk = 2 * min(sum(permuted_statistics <= twelvewk_sus_rem_observedratio), sum(permuted_statistics >= twelvewk_sus_rem_observedratio)) / num_permutations;
end

fprintf('Two-tailed p-value: %.4f\n', p_value_twowk);
fprintf('Two-tailed p-value: %.4f\n', p_value_fourwk);
fprintf('Two-tailed p-value: %.4f\n', p_value_sixwk);
fprintf('Two-tailed p-value: %.4f\n', p_value_eightwk);
fprintf('Two-tailed p-value: %.4f\n', p_value_twelvewk);

figure
histogram(permuted_statistics, 'NumBins', 150);
