close all
clear all

%Import file 
optsA = detectImportOptions('Three-D_HAMDTotals.xlsx');
optsA.VariableNames = {'StudyID','HAMD-0','HAMD-1','HAMD-2','HAMD-3','HAMD-4','HAMD-5','HAMD-End','HAMD-1F','HAMD-2F','HAMD-3F'};
optsA.VariableTypes = {'string','double','double','double','double','double','double','double','double','double','double'};
HAMD_totals = readtable('Three-D_HAMDTotals.xlsx',optsA);

% %Replace NaN with zeros
% numericData = HAMD_totals{:, 2:end};
% nanIndices = isnan(numericData);
% numericData(nanIndices) = 0;
% HAMD_totals(:, 2:end) = array2table(numericData);

HAMD_totals_array = table2array(HAMD_totals);
%HAMD_totals_array = str2double(HAMD_totals_array);
HAMD_totals_array_cropped = HAMD_totals_array(:,2:end);
HAMD_totals_array_cropped = str2double(HAMD_totals_array_cropped);
%% Initialize variables

meetsCriteria = false(size(HAMD_totals_array_cropped)); 
allHAMD_Delta = zeros(size(HAMD_totals_array_cropped)); 
sg_reversalvalue = [];

%% Evaluate sudden gains for each ppt
% Loop through each participant
for sessionindex = 2:8
    % Compute differences between consecutive scores for each participant
    HAMD_Delta =  HAMD_totals_array_cropped(:, sessionindex + 1) - HAMD_totals_array_cropped(:, sessionindex);

    allHAMD_Delta(:, sessionindex) = HAMD_Delta;

    halfdelta = allHAMD_Delta./2; %divide the change in HAMD by 2

    reversalvalue = HAMD_totals_array_cropped + halfdelta;

    % Calculate the percentage differences
    HAMD_percentdecrease =  (HAMD_totals_array_cropped(:, sessionindex) - HAMD_totals_array_cropped(:, sessionindex + 1)) ./ HAMD_totals_array_cropped(:, sessionindex) * 100;

    % Extract 3 Pre and 3 Post HAMD values from gain
    if sessionindex == 2
         HAMD_pre = HAMD_totals_array_cropped(:, (sessionindex - 1):(sessionindex)); 
    else 
        HAMD_pre = HAMD_totals_array_cropped(:, (sessionindex - 2):(sessionindex)); 
    end

    if sessionindex == 8
         HAMD_post = HAMD_totals_array_cropped(:, (sessionindex + 1):(sessionindex + 2)); 
    else 
        HAMD_post = HAMD_totals_array_cropped(:, (sessionindex + 1):(sessionindex + 3));
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
            symptomfluctuation = 4.303 * sqrt(((HAMD_totals_array_cropped(PrePost_index, sessionindex) - 1) .* std_pre(PrePost_index).^2 + (HAMD_totals_array_cropped(PrePost_index, sessionindex + 1) - 1) .* std_post(PrePost_index).^2) ./ (HAMD_totals_array_cropped(PrePost_index, sessionindex) + HAMD_totals_array_cropped(PrePost_index, sessionindex + 1) - 2));

            allsymptomfluctuation(PrePost_index,:) = symptomfluctuation;
        end 

        % Apply symptom fluctuation formula for one missing value 
        for j = 1:length(One_NaN_rows)
            OneNaN_index = One_NaN_rows(j);
            symptomfluctuation = 3.182 * sqrt(((HAMD_totals_array_cropped(OneNaN_index, sessionindex) - 1) .* std_pre(OneNaN_index).^2 + (HAMD_totals_array_cropped(OneNaN_index, sessionindex + 1) - 1) .* std_post(OneNaN_index).^2) ./ (HAMD_totals_array_cropped(OneNaN_index, sessionindex) + HAMD_totals_array_cropped(OneNaN_index, sessionindex + 1) - 2));

            allsymptomfluctuation(OneNaN_index,:) = symptomfluctuation;
        end 

        % Apply symptom fluctuation formula for no missing values
        for k = 1:length(zero_rows);
            zero_index = zero_rows(k);
        symptomfluctuation = 2.776 * sqrt(((HAMD_totals_array_cropped(zero_index, sessionindex) - 1) .* std_pre(zero_index).^2 + (HAMD_totals_array_cropped(zero_index, sessionindex + 1) - 1) .* std_post(zero_index).^2) ./ (HAMD_totals_array_cropped(zero_index, sessionindex) + HAMD_totals_array_cropped(zero_index, sessionindex + 1) - 2));

            allsymptomfluctuation(zero_index,:) = symptomfluctuation;
        end

        % Criteria 3: Calculate Mpre - Mpost
    MeanDelta = mean(HAMD_pre, 2, 'omitmissing') - mean(HAMD_post, 2, 'omitmissing');
    
    %Check if all three criteria are met for each score difference
    criteriaCheck = HAMD_Delta <= -7 & HAMD_percentdecrease >= 25 &  MeanDelta > allsymptomfluctuation;

    meetsCriteria(:, sessionindex) = criteriaCheck;
  
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
consolidatedreversal = false(size(HAMD_totals_array_cropped));

% Compare HAMD values collected after the sudden gain to the reversal value
% rowIndex indicates rows where sudden gain occured
for i = 1:numel(rowIndex);

    % Store indices of values after the sudden gain
    indices = zero_after_one_indices{rowIndex(i)};
    
    % sg_index(i,2) = column index where sudden gain even toccured
    indices = indices(indices >= sg_index(i,2));

    % If the HAMD is greater than reversal value, store in reversal_ppt
reversal_ppt = HAMD_totals_array_cropped(rowIndex(i), indices) >= sg_index(i,3);

consolidatedreversal(rowIndex(i), indices) = reversal_ppt;
end 

% Identify rows that meet all criteria
rows_with_ones_sg = any(meetsCriteria == 1, 2);
rows_with_ones_reversals = any(consolidatedreversal == 1, 2);

indices_with_ones = find(rows_with_ones_sg);
indices_with_ones_reversals = find(rows_with_ones_reversals);

% Create a new variable and assign modified values
rows_with_ones_sg_reversed = rows_with_ones_sg;
rows_with_ones_sg_reversed(indices_with_ones_reversals) = 0;

% %% Find sudden gain participants and which session they experience their first sudden gain event
% 
% final_non_sgppt = find(rows_with_ones_sg_reversed == 0);
% 
% % Create matrix that only keeps the first sudden gain and deletes second
% % sudden gain event for the same participant
% 
% % Initialize the modified matrix with zeros
% modifiedmeetsCriteria = zeros(size(meetsCriteria));
% 
% % Loop through each row
% for i = 1:size(meetsCriteria, 1)
%     % Find indices of 1's in the current row
%     indicesOfOnes = find(meetsCriteria(i, :) == 1);
% 
%     % Keep the first 1 (if any)
%     if ~isempty(indicesOfOnes)
%         modifiedmeetsCriteria(i, indicesOfOnes(1)) = 1;
% 
%         % Make the second 1 (if any) a 0
%         if numel(indicesOfOnes) > 1
%             modifiedmeetsCriteria(i, indicesOfOnes(2:end)) = 0;
%         end
%     end 
% 
% end
% 
%     modifiedmeetsCriteria(final_non_sgppt, :) = 0;
% 
%       %Number of patients who experienced their first sudden gain at each
%     %session
%     sgppt_across_sessions = sum(modifiedmeetsCriteria, 1);
% 
%     %Sums up all the people who have experienced a sudden gain prior to
%     %each session
%     cumulativesum_sgppt = cumsum(sgppt_across_sessions);
% 
% 
% sgavg_across_sessions = sum(meetsCriteria, 1); %sum total sudden gains across each session for every participant

 % %Number of patients who experienced their first sudden gain at each
 %    %session
 %    sgppt_across_sessions = sum(modifiedmeetsCriteria, 1);
 % 
 %    %Sums up all the people who have experienced a sudden gain prior to
 %    %each session
 %    cumulativesum_sgppt = cumsum(sgppt_across_sessions);
 % 
 %    shifted_cumulativesum_sgppt = [0, cumulativesum_sgppt(1:18)];


 %% Initialize the number of rows and columns for 4x4 subplotx
numRows = 5;
numCols = 6;

% Calculate the total number of subplots
numSubplots = numRows * numCols;

% Create a cell array to hold the handles to individual figures
figureHandles = cell(1, numSubplots);

%% Generate figures for each participant

variableNames = {'sg_ppt', 'remission_met_ppt'};

% Convert indicies_with_ones from a logical to a table format
rows_with_ones_table = array2table(rows_with_ones_sg_reversed, 'variableNames', {'sg_ppt'});

% Concatenate sg_ppt to HAMD_by_session
HAMD_totals = [HAMD_totals, rows_with_ones_table];

% Get the number of people in the table
TotalParticipants = size(HAMD_totals_array_cropped, 1);

% Sudden Gain Analysis Comparison: Nili vs Us
our_sgsum = sum(HAMD_totals.sg_ppt == 1);
allour_gainers = (HAMD_totals.sg_ppt == 1);

% Run some preliminary statistics on SGs and compare to Nili's findings
our_sg_proportion = (our_sgsum / 240) * 100;
fprintf('Our total sudden gainers: %d\n', our_sgsum);
%fprintf('Percentage of participants who experienced a sudden gain: %.2f%%\n', our_sg_proportion)


% Loop through each person and create a figure
for personIndex = 1:TotalParticipants
    % Store all HAMD scores for the current person into an array
    individual_HAMD = HAMD_by_session(personIndex, 2:19);
    sg_ppt = HAMD_by_session.sg_ppt(personIndex);

    %% NOT EDITED YET
    % Determine the figure and subplot index
    figureIndex = ceil(personIndex / numSubplots);
    subplotIndex = mod(personIndex - 1, numSubplots) + 1;

    % If it's the first subplot of a new figure, create a new figure
    if subplotIndex == 1
        figureHandles{figureIndex} = figure;
    end

    % Use subplot to specify the current subplot within the current figure
    subplot(numRows, numCols, subplotIndex);

    % % Use subplot to specify the current subplot
    % subplot(numRows, numCols, mod(personIndex - 1, numSubplots) + 1);

    %Index through the ID column
    ID = HAMD_by_session.ID(personIndex);

      %Plot the HAMD data for the current person
    if sg_ppt == 1 && remission_met_ppt == 1 % all_crit == 1
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-m*')

        % Add labels, titles, legends, etc. to the figure
        title(['Subject ', num2str(ID), ' SG: session ', num2str(sg_session_n), ' R: ', num2str(remission_session)]);
        xlabel('ECT Treatment Visit');
        ylabel('HAM-D Score');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(60)]);

%         hold on;
% 
% % Plotting specific points in a different color based on the logical array
% plot(0:(width(individual_HAMD) - 1), HAMD_by_session(meetsCriteria), 'ro', 'MarkerSize', 8); % Plot red circles for specific points
% 
% hold off;

        
    elseif sg_ppt == 1 && remission_met_ppt == 0 % all_crit == 0 
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-r*')

        % Add labels, titles, legends, etc. to the figure
        title(['Subject ', num2str(ID), ' SG: session ', num2str(sg_session_n)]);
        xlabel('ECT Treatment Visit');
        ylabel('HAM-D Score');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(60)]);

%         hold on;
% 
% % Plotting specific points in a different color based on the logical array
% plot(0:(width(individual_HAMD) - 1), HAMD_by_session(meetsCriteria), 'ro', 'MarkerSize', 8); % Plot red circles for specific points
% 
% hold off;
        
    elseif sg_ppt == 0 && remission_met_ppt == 1 % all_crit == 1
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-b*')

        % Add labels, titles, legends, etc. to the figure
        title(['Subject ', num2str(ID), ' R: ', num2str(remission_session)]);
        xlabel('ECT Treatment Visit');
        ylabel('HAM-D Score');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(60)]);

%         hold on;
% 
% % Plotting specific points in a different color based on the logical array
% plot(0:(width(individual_HAMD) - 1), HAMD_by_session(meetsCriteria), 'ro', 'MarkerSize', 8); % Plot red circles for specific points
% 
% hold off;
        
    else plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-k*')

    % Add labels, titles, legends, etc. to the figure
    title(['Subject ', num2str(ID)]);
    xlabel('ECT Treatment Visit');
    ylabel('HAM-D Score');

    %Set xlim and ylim
    xlim([0, max(18)]);
    ylim([min(0), max(60)]);

%     hold on;
% 
% % Plotting specific points in a different color based on the logical array
% plot(0:(width(individual_HAMD) - 1), HAMD_by_session(meetsCriteria), 'ro', 'MarkerSize', 8); % Plot red circles for specific points
% 
% hold off;

    end
end

%% Remission & Sudden Gain Statistics & Comparison
% Extract the instances of remission and sudden gain for each ppt
remission_met_1 = HAMD_by_session.remission_met_ppt(HAMD_by_session.remission_met_ppt ~= 0);
sg_ppt_1 = HAMD_by_session.sg_ppt(HAMD_by_session.sg_ppt ~= 0);

HAMD_by_session_array_all = table2array(HAMD_by_session);

% Find ppt who were both a remitter and sudden gainer
rem_and_sg_instances = all(HAMD_by_session.sg_ppt == 1 & HAMD_by_session.remission_met_ppt == 1, 2);

% Calculate total of ppts who were both remitters and sudden gainers, find percent of total participants
sum_rem_and_sg = sum(rem_and_sg_instances);
rem_sg_percentage = (sum_rem_and_sg / 240) * 100;
rem_and_sg_indices = find(rem_and_sg_instances == 1);

% Calculate proportion of remitters that were also sudden gainers
sg_of_remitters = (sum(rem_and_sg_instances) / numel(remission_met_1)) * 100;

% Calculate proportion of sudden gainers that were also remitters
remitters_of_sg = (sum(rem_and_sg_instances) / numel(sg_ppt_1)) * 100;

% Display all previously calculated statistics
fprintf('%d participants were both remitters and sudden gainers\n', sum_rem_and_sg)
fprintf('%.2f%% of participants were both remitters and sudden gainers\n', rem_sg_percentage)
fprintf('%.2f%% of remitters were sudden gainers\n', sg_of_remitters)
fprintf('%.2f%% of sudden gainers were remitters\n', remitters_of_sg)

%% Find denominators for normalized SG histogram
% Initialize variable
num_remaining_ppt = zeros(size(sgavg_across_sessions));

% Loop through each session 1-18 to find # of remaining ppts at each session
for sessionIndex = 2:19
    num_remaining_ppt(sessionIndex) = sum(~isnan(HAMD_by_session_array(:, sessionIndex)));
end

%% Plot proportion of people who experienced sudden gains at each session over remaining total number of participants
       figure
       sgavg_over_remainingppt = sgavg_across_sessions./num_remaining_ppt;
   bar(1:(width(individual_HAMD)), sgppt_over_remainingppt(2:19), 'b');


        % Add labels, titles, legends, etc. to the figure
        title(['Sudden Gain Events Each Session Over Total Remaining Participants']);
        xlabel('ECT Treatment Visit');
        ylabel('Proportion of Participants Over Total Remaining Participants');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(0.15)]);

%% Plot proportion of people who experienced sudden gains at each session over remaining total number of participants
       figure
       sgppt_over_remainingppt = sgppt_across_sessions./num_remaining_ppt;
       % sgppt_over_remainingppt = sgavg_across_sessions(1:19)./num_remaining_ppt;
   bar(1:(width(individual_HAMD)), sgppt_over_remainingppt(2:19), 'b');


        % Add labels, titles, legends, etc. to the figure
        title(['Sudden Gain Participants Who Experience Their First Sudden Gain Event At Each Session Over Total Remaining Participants']);
        xlabel('ECT Treatment Visit');
        ylabel('Proportion of Participants Over Total Remaining Participants');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(0.15)]);

%% Plot quantity of sudden gains across sessions
figure
   %plot(0:(width(individual_HAMD)), sgavg_across_sessions(1:19), '*')
   bar(1:(width(individual_HAMD)), sgavg_across_sessions(2:end), 'b');

        % Add labels, titles, legends, etc. to the figure
        title(['Sudden Gain Events Across Sessions']);
        xlabel('ECT Treatment Visit');
        ylabel('Number of Sudden Gains');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(32)]);

