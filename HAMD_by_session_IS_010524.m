% % Transform Nili's Outcome and Static Predictor Table to Continuoius
% Process form % % 
clear all
%% load table
opts = delimitedTextImportOptions("NumVariables", 88);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ID", "H_TOTALSCORE_0", "H_TOTALSCORE_1", "H_TOTALSCORE_2", "H_TOTALSCORE_3", "H_TOTALSCORE_4", "H_TOTALSCORE_5", "H_TOTALSCORE_6", "H_TOTALSCORE_7", "H_TOTALSCORE_8", "H_TOTALSCORE_9", "H_TOTALSCORE_10", "H_TOTALSCORE_11", "H_TOTALSCORE_12", "H_TOTALSCORE_13", "H_TOTALSCORE_14", "H_TOTALSCORE_15", "H_TOTALSCORE_16", "H_TOTALSCORE_17", "sg_session_n", "all_crit", "remission_met", "remission_session"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";

% Import the data
HAMD_by_session = readtable("/Users/ins4004/Documents/Projects/ECT_Response_Trajectories_Project/PRIDE_Data_Analyses/Derived/id_HAMDtotal", opts);
HAMD_by_session_array = table2array(HAMD_by_session);
HAMD_by_session_array(:, 20:23) = [];

%% Detect sudden gains

% Initialize a logical array to store results
meetsCriteria = false(size(HAMD_by_session_array)); % Assuming 240 rows
 allHAMD_Delta = zeros(size(HAMD_by_session_array)); 
 sg_reversalvalue = [];

% % Get the number of people in the table
% TotalSessions = 2:width(HAMD_by_session);

% Initialize logical array to store remission data starting at session 2
remission_criteria_met = false(size(HAMD_by_session_array, 1), 16);

% Loop through each participant and extract their data
for p = 1:size(HAMD_by_session_array, 1)
    participant_data = HAMD_by_session_array(p, :);
    % Create value to check if two consecutive HAMD values are <= 10, and that any consecutive increase is either <= 3, or value of <= 6
    remission_criteria_met_10 = (participant_data(3:end-1) <= 10) & (participant_data(4:end) <= 10);
    increase = participant_data(4:end) - participant_data(3:end-1);
    second_less_than_6 = participant_data(4:end) <= 6;
    remission_criteria_met_ppt = remission_criteria_met_10 & ((increase <= 3) | second_less_than_6);
    % Check all remission criteria for each participant
    remission_criteria_met(p, :) = remission_criteria_met_ppt;
end

%% Check if remission criteria is met for each ppt
% Initialize variable to store remission criteria data
remission_criteria_met_atall = false(size(remission_criteria_met, 1), 1);

% Check if remission criteria is met at least once for each ppt
remission_criteria_met_atall = any(remission_criteria_met == 1, 2);

% Get sum of total remitters
total_remitters = sum(remission_criteria_met_atall);

%% Index each participants first instance of remission
% Extract the row (ppt) and column (session) that remission was reached for each instance of remission
[row_remission_index, col_remission_index] = find(remission_criteria_met == 1);

% Correct the session data to align with the original data by adding 1 the indexed sessions
col_remission_index_corrected = col_remission_index + 1;

% Concatenize the ppt and session data into the same 2-column variable and sort it so that it is in order through the ppts
remission_instances = horzcat(row_remission_index, col_remission_index_corrected);
remission_instances_sorted = sortrows(remission_instances, 1);

% Extract the first instances of remission for each ppt
[first_rem_instance, rem_sessionidx] = unique(remission_instances_sorted(:, 1), 'stable');
first_rem_session = remission_instances_sorted(rem_sessionidx, 2);

% Concatenize first remission instance data into one 2-column variable
remission_met_indices = horzcat(first_rem_instance, first_rem_session);

%% Evaluate sudden gains for each ppt
% Loop through each participant
for sessionindex = 3:17
    % Compute differences between consecutive scores for each participant
    HAMD_Delta =  HAMD_by_session_array(:, sessionindex + 1) - HAMD_by_session_array(:, sessionindex);

    allHAMD_Delta(:, sessionindex) = HAMD_Delta;

    halfdelta = allHAMD_Delta./2; %divide the change in HAMD by 2

    reversalvalue = HAMD_by_session_array + halfdelta;

    % Calculate the percentage differences
    HAMD_percentdecrease =  (HAMD_by_session_array(:, sessionindex) - HAMD_by_session_array(:, sessionindex + 1)) ./ HAMD_by_session_array(:, sessionindex) * 100;

    % Extract 3 Pre and 3 Post HAMD values from gain
    if sessionindex == 3
         HAMD_pre = HAMD_by_session_array(:, (sessionindex - 1):(sessionindex)); 
    else 
        HAMD_pre = HAMD_by_session_array(:, (sessionindex - 2):(sessionindex)); 
    end

    if sessionindex == 17
         HAMD_post = HAMD_by_session_array(:, (sessionindex + 1):(sessionindex + 2)); 
    else 
        HAMD_post = HAMD_by_session_array(:, (sessionindex + 1):(sessionindex + 3));
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
  
    % % Check for overlapped indexes of participants
    % commonIndexes = intersect(One_NaN_rows, PrePost_NaN_rows);
    % commonIndexes2 = intersect(PrePost_NaN_rows, bye_rows);

    % Create a range of all possible row indexes
    totalRows = 1:240; % Replace 'totalNumberOfRows' with the total number of rows in your data

    % Combine all indexes into one array
    allIndexes = sort([zero_rows; One_NaN_rows; PrePost_NaN_rows; bye_rows]);

    % Find the missing rows
    missingRows = setdiff(totalRows, unique(allIndexes));

    doublecheck = length(One_NaN_rows) + length(PrePost_NaN_rows) + length(bye_rows) + length(zero_rows);

    % Apply criteria 3, adjusting for missing values
            
    allsymptomfluctuation(bye_rows,:) = NaN;
       

        for i = 1:length(PrePost_NaN_rows)
            PrePost_index = PrePost_NaN_rows(i);
            symptomfluctuation = 4.303 * sqrt(((HAMD_by_session_array(PrePost_index, sessionindex) - 1) .* std_pre(PrePost_index).^2 + (HAMD_by_session_array(PrePost_index, sessionindex + 1) - 1) .* std_post(PrePost_index).^2) ./ (HAMD_by_session_array(PrePost_index, sessionindex) + HAMD_by_session_array(PrePost_index, sessionindex + 1) - 2));

            allsymptomfluctuation(PrePost_index,:) = symptomfluctuation;
        end 

        for j = 1:length(One_NaN_rows)
            OneNaN_index = One_NaN_rows(j);
            symptomfluctuation = 3.182 * sqrt(((HAMD_by_session_array(OneNaN_index, sessionindex) - 1) .* std_pre(OneNaN_index).^2 + (HAMD_by_session_array(OneNaN_index, sessionindex + 1) - 1) .* std_post(OneNaN_index).^2) ./ (HAMD_by_session_array(OneNaN_index, sessionindex) + HAMD_by_session_array(OneNaN_index, sessionindex + 1) - 2));

            allsymptomfluctuation(OneNaN_index,:) = symptomfluctuation;
        end 

        for k = 1:length(zero_rows);
            zero_index = zero_rows(k);
        symptomfluctuation = 2.776 * sqrt(((HAMD_by_session_array(zero_index, sessionindex) - 1) .* std_pre(zero_index).^2 + (HAMD_by_session_array(zero_index, sessionindex + 1) - 1) .* std_post(zero_index).^2) ./ (HAMD_by_session_array(zero_index, sessionindex) + HAMD_by_session_array(zero_index, sessionindex + 1) - 2));

        allsymptomfluctuation(zero_index,:) = symptomfluctuation;
        end

        % Criteria 3: Calculate Mpre - Mpost
    MeanDelta = mean(HAMD_pre, 2, 'omitmissing') - mean(HAMD_post, 2, 'omitmissing');
    
    %Check if all three criteria are met for each score difference
    criteriaCheck = HAMD_Delta <= -7 & HAMD_percentdecrease >= 25 &  MeanDelta > allsymptomfluctuation;

    meetsCriteria(:, sessionindex) = criteriaCheck;
  
end

%% Find reversals occuring for each participant
[rowIndex, colIndex] = find(meetsCriteria == 1);

sg_reversalvalue = reversalvalue(sub2ind(size(meetsCriteria), rowIndex, colIndex));

sg_index = horzcat(rowIndex, colIndex, sg_reversalvalue);

sg_linearindex = sub2ind(size(meetsCriteria), rowIndex, colIndex);

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
    zero_after_one_indices{row} = unique(zero_indices);
end


% % Indices after sudden gains: Convert row/column indices into linear indices
% for i = 1:numel(zero_after_one_indices);
%     % Process each cell and convert subscripts to linear indices
%     % Replace 'matrix_size' with the actual size of your matrix
%     matrix_size = [240, 19]; % Replace nRows and nCols with actual dimensions
%     incidents_aftersg_linearindex{i} = sub2ind(matrix_size, repmat(i, size(zero_after_one_indices{i})), zero_after_one_indices{i});
% end
% 
% % Initialize an empty array to store consolidated values
% indexconsolidated_aftersg = [];
% 
% % Loop through each cell
% for i = 1:numel(incidents_aftersg_linearindex)
%     % Check if the cell is not empty and is an array
%     if ~isempty(incidents_aftersg_linearindex{i}) 
%         % Concatenate the array in the cell to the consolidated array
%         indexconsolidated_aftersg = [indexconsolidated_aftersg, incidents_aftersg_linearindex{i}];
%     end
% end

% indexconsolidated_aftersg = indexconsolidated_aftersg';

consolidatedreversal = false(size(HAMD_by_session_array));

% Indicate reversal incident that occurs by session and participant
for i = 1:numel(rowIndex);

    indices = zero_after_one_indices{rowIndex(i)};
    
    % Exclude values under the threshold
    indices = indices(indices >= sg_index(i,2));

reversal_ppt = HAMD_by_session_array(rowIndex(i), indices) >= sg_index(i,3);

consolidatedreversal(rowIndex(i), indices) = reversal_ppt;
end 


sgavg_across_sessions = sum(meetsCriteria, 1); %sum total sudden gains across each session for every participant

% Identify rows that meet all criteria
rows_with_ones_sg = any(meetsCriteria == 1, 2);
rows_with_ones_reversals = any(consolidatedreversal == 1, 2);

indices_with_ones = find(rows_with_ones_sg);
indices_with_ones_reversals = find(rows_with_ones_reversals);

% Create a new variable and assign modified values
rows_with_ones_sg_reversed = rows_with_ones_sg;
rows_with_ones_sg_reversed(indices_with_ones_reversals) = 0;


 %% Initialize the number of rows and columns for 4x4 subplotx
numRows = 4;
numCols = 6;

% Calculate the total number of subplots
numSubplots = numRows * numCols;

% Create a cell array to hold the handles to individual figures
figureHandles = cell(1, numSubplots);

%% Generate figures for each participant

variableNames = {'sg_ppt', 'remission_met_ppt'};

% Convert indicies_with_ones from a logical to a table format
rows_with_ones_table = array2table(rows_with_ones_sg_reversed, 'variableNames', {'sg_ppt'});


% Convert remission_criteria_met_atall from a logical to a table format
remission_criteria_met_atall_table = array2table(remission_criteria_met_atall, 'variableNames', {'remission_met_ppt'});

% Concatenate sg_ppt to HAMD_by_session
HAMD_by_session = [HAMD_by_session, rows_with_ones_table];

HAMD_by_session = [HAMD_by_session, remission_criteria_met_atall_table];

% Get the number of people in the table
TotalParticipants = size(HAMD_by_session, 1);

% Sudden Gain Analysis Comparison: Nili vs Us
our_sgsum = sum(HAMD_by_session.sg_ppt == 1)
Nili_sgsum = sum(HAMD_by_session.all_crit == 1)
allour_gainers = (HAMD_by_session.sg_ppt == 1);
matching_gainers = all(HAMD_by_session.sg_ppt == 1 & HAMD_by_session.all_crit == 1, 2);

onlyour_gainers = all(HAMD_by_session.sg_ppt == 1 & HAMD_by_session.all_crit == 0, 2);
onlyNili_gainers = all(HAMD_by_session.sg_ppt == 0 & HAMD_by_session.all_crit == 1, 2);
sum_matching_gainers = sum(matching_gainers == 1)

row_onlyour_gainers = find(onlyour_gainers);
row_onlyNiligainers = find(onlyNili_gainers);
row_matchinggainers = find(matching_gainers);
row_ourgainers = find(allour_gainers);

% Run some preliminary statistics on SGs and compare to Nili's findings
our_sg_proportion = (our_sgsum / 240) * 100;
fprintf('Our total sudden gainers: %d\n', our_sgsum);
fprintf('Nili total sudden gainers: %d\n', Nili_sgsum);
fprintf('Total matching sudden gainers: %d\n', sum_matching_gainers);
fprintf('Percentage of participants who experienced a sudden gain: %.2f%%\n', our_sg_proportion)


%% Run statistics on remission
% Find amount of people who met remission and which sessions they met remission based on all.7 sheet
sum_remitted = sum(HAMD_by_session.remission_met_ppt == 1);

% Calculate the remission rate
remission_rate = (sum_remitted / 240) * 100;

% Find most common session for remission to occur

% Extract nonremitters from session count
nonZerosessions = HAMD_by_session.remission_session(HAMD_by_session.remission_session ~= 0);

% Find most common session
most_common_remission_session = mode(nonZerosessions);
unique_sessions = unique(nonZerosessions);

% Count amount of remission instances at each session
counts = histcounts(nonZerosessions);
counts_vert = counts';

% Combine session indices with counts
remission_session_freq = horzcat(unique_sessions, counts_vert);

% Find percentage of remission instances that occur at each session
remission_session_percentage = (counts_vert ./ 144) * 100;

% Concatenize session, # of remitters, and percentage of remitters into one 144x3 variable
remission_session_freq = horzcat(unique_sessions, counts_vert, remission_session_percentage);

% if remission_met == 1
%     sg_of_remitters = find(sg_ppt == 1);
%     fprintf('%d of remitters were sudden gainers', sg_of_remitters)
% end

% Prints the number, percentage, session mode, and breakdown of frequency of remission by session
fprintf('Total remitters: %d\n', sum_remitted);
fprintf('Remission rate: %d%%\n', remission_rate);
fprintf('Most remitted by session: %d\n', most_common_remission_session);
for m = 1:length(unique_sessions)
    fprintf('Session: %d: Frequency %.2f%%\n', remission_session_freq(m, 1), remission_session_freq(m, 3))
end


% Loop through each person and create a figure
for personIndex = 1:TotalParticipants
    % Store all HAMD scores for the current person into an array
    individual_HAMD = HAMD_by_session(personIndex, 2:19);
    all_crit = HAMD_by_session.all_crit(personIndex);
    sg_session_n = HAMD_by_session.sg_session_n(personIndex);
    sg_ppt = HAMD_by_session.sg_ppt(personIndex);
    remission_met_ppt = HAMD_by_session.remission_met_ppt(personIndex);
    remission_session = HAMD_by_session.remission_session(personIndex);
    
    % if remission_met && sg_ppt == 1
    %     rem_sg(personIndex) = true;
    % else
    %     rem_sg(personIndex) = false;
    % end
    % rem_sg_vert = rem_sg';
    % sum_rem_sg_vert = sum(rem_sg_vert);
    % rem_sg_percentage = (sum_rem_sg_vert / 240) * 100; 
    % sg_of_remitters = sum(sg_ppt) / sum(remission_met);

    % fprintf('%.2f%%\n of participants were both remitters and sudden gainers', rem_sg_percentage)

    % % Create a figure for the current person
    % figure;

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

figure

%% Plot quantity of sudden gains across sessions
   %plot(0:(width(individual_HAMD)), sgavg_across_sessions(1:19), '*')
   bar(1:(width(individual_HAMD)), sgavg_across_sessions(2:end), 'b');

        % Add labels, titles, legends, etc. to the figure
        title(['Sudden Gains Across Sessions']);
        xlabel('ECT Treatment Visit');
        ylabel('Number of Sudden Gains');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(32)]);

%     % If we have created a full set of subplots, create a new figure
%     if mod(personIndex, numSubplots) == 0
%         for i = 1:numSubplots
%             subplot(numRows, numCols, i);
%             figureHandles(i) = get(gca, 'Children');
%         end
%         % saveas(gcf, ['SubplotFigure', num2str(personIndex / numSubplots), '.png']);  % Save the subplot as a single figure
%         % close all;  % Close all figures to start a new subplot set
%     end
% end
% 
% % Save the last subplot set as a single figure
% for i = 1:numSubplots
%     subplot(numRows, numCols, i);
%     figureHandles(i) = get(gca, 'Children');
% end
% % saveas(gcf, ['SubplotFigure', num2str(personIndex / numSubplots), '.png']);
% % close all;  % Close all figures

%% Detect any instances of HAMD <= 10
% remission_crit1 = zeros(size(HAMD_by_session_array));
% [rowHAMD_under10, colHAMD_under10] = find(HAMD_by_session_array <= 10);
% remission_crit1_index = horzcat(rowHAMD_under10, colHAMD_under10);
% remission_crit1_indexsorted = sortrows(remission_crit1_index);
% [first_instance, sessionidx] = unique(remission_crit1_indexsorted(:, 1), 'stable');
% first_instance_session = remission_crit1_indexsorted(sessionidx, 2);
% remission_crit1_met_indices = horzcat(first_instance, first_instance_session);
%% Store the next value after first session where HAMD <= 10
% Initialize overall remission criteria 1 variable
% remission_crit1_overall = [];
% remission_crit1_post10_indices = zeros(length(remission_crit1_met_indices), 1);
% remission_crit1_overall_indices = zeros(length(remission_crit1_met_indices), 3);

% Loop through each participant with a HAMD <= 10 and store the indices of the next consecutive score
% for i = 1:length(remission_crit1_met_indices)
% 
%     remission_crit1_post10_indices(i) = remission_crit1_met_indices(i, 2) + 1;
% 
%     
% 
%     remission_crit1_overall_indices = horzcat(first_instance, first_instance_session, remission_crit1_post10_indices);
% end

% Initialize variable for storing the HAMD value at the next consecutive score after HAMD <= 10
% remission_crit1_second_value = zeros(size(remission_crit1_overall_indices, 1), 1);
% remission_crit1_overall_met = false(size(remission_crit1_overall_indices, 1), 1);
% remission_crit1_notmet_check_second_indices = zeros(size(remission_crit1_overall_indices, 1), 1);
% remission_crit1_all_met = false(size(HAMD_by_session_array, 1), 1);


% Loop through each participant with a HAMD <= 10 and store the next consecutive HAMD score
% for j = 1:length(remission_crit1_overall_indices)
% 
%     remission_crit1_second_value(j) = HAMD_by_session_array(remission_crit1_overall_indices(j, 1), remission_crit1_overall_indices(j, 3));
    
    % remission_crit1_data = horzcat(first_instance, first_instance_session, remission_crit1_post10_indices, remission_crit1_second_value);
    
    % Check if next consecutive HAMD score is <= 10
    % remission_crit1_overall_met(j) = remission_crit1_second_value(j) <= 10;
    % 
    % remission_crit1_data = horzcat(first_instance, first_instance_session, remission_crit1_post10_indices, remission_crit1_second_value, remission_crit1_overall_met);
    % 
    % if remission_crit1_data(j, 5) == 0
    %     remission_crit1_notmet_check_second_indices(j) = remission_crit1_overall_indices(j, 3) + 1;
    % 
    %     remission_crit1_data = horzcat(first_instance, first_instance_session, remission_crit1_post10_indices, remission_crit1_second_value, remission_crit1_overall_met, remission_crit1_notmet_check_second_indices);
    % 
    %    remission_crit1_notmet_check_second_indices_nonzero = remission_crit1_data(:, 6) ~= 0;
    % 
    %    remission_crit1_notmet_check_second_indices_nonzero_values = remission_crit1_data(remission_crit1_notmet_check_second_indices_nonzero, 6);
    % 
    %    remission_crit1_notmet_check_second_indices_nonzero_ppt = remission_crit1_data(remission_crit1_notmet_check_second_indices_nonzero, 1);
    % 
    %    remission_crit1_second_indices_ppt = horzcat(remission_crit1_notmet_check_second_indices_nonzero_ppt, remission_crit1_notmet_check_second_indices_nonzero_values);
    % 
    %     num_elements = size(remission_crit1_notmet_check_second_indices_nonzero_ppt, 1);
    % 
    %     HAMD_values_post_10 = cell(num_elements, 1);
    % 
    %     for k = 1:num_elements
    %          HAMD_values_post_10{k} = HAMD_by_session_array(remission_crit1_second_indices_ppt(k, 1), remission_crit1_notmet_check_second_indices_nonzero_values(k)+1:end);
    %          less_than_10 = cellfun(@(x) any(x <= 10 & [x(2:end) <= 10, false]), HAMD_values_post_10);
    %          remission_crit1_second_indices_ppt = horzcat(remission_crit1_notmet_check_second_indices_nonzero_ppt, remission_crit1_notmet_check_second_indices_nonzero_values, less_than_10);
    %     end
       
        % current_index = remission_crit1_data(j, 1);
        % 
        % % remission_crit1_all_met(remission_crit1_data(j, 1)) = remission_crit1_overall_met(remission_crit1_data(j, 1));
        % 
        % % Check if the current index exists in the smaller set
        % idx_43 = find(remission_crit1_second_indices_ppt(:, 1) == current_index);
        % 
        % if ~isempty(idx_43)
        %     % If the index is found in the smaller set, assign its value
        %     remission_crit1_all_met(j) = less_than_10(idx_43);
        % else
        %     % If the index is not found in the smaller set, assign value from the larger set
        %     remission_crit1_all_met(j) = remission_crit1_overall_met(j);
        % end
        % 

%     end
% end

% remission_crit1_second_value_vert = remission_crit1_second_value';
% remission_crit1_post10_indices = remission_crit1_post10_indices';

%remission_crit1_overall = HAMD_by_session_array(remission_crit1_post10_indices(i)) <= 10;




% remission_crit1 = all(HAMD_by_session_array, 2:end) <= 10)

% initialize first instant indices
% first_instance_under10 = zeros(240, 1);
% 
% % Get the number of people in the table
% TotalParticipants = size(HAMD_by_session_array, 1);
% 
% % Detect first instances of HAMD <= 10 for each participant
% for personindex = 1:TotalParticipants
%     indices_less_than_10 = find(HAMD_by_session_array(personindex, :) <= 10);
%     if ~isempty(indices_less_than_10)
%         first_instance_under10(personindex) = indices_less_than_10(1);
%     else
%         first_instance_under10(personindex) = 0;
%     end
% 
% end
% 
% rows_with_remission = first_instance_under10 ~= 0;
% for personindex = 1:TotalParticipants
%     if rows_with_remission(personindex) == 1
%         remission_crit1_sessionindex = HAMD_by_session_array(:, first_instance_under10(personindex));
%     else
%         remission_crit1_sessionindex = 0;
%     end
% end



% remission_crit2_sessionindex = first_instance_under10(personindex, remission_crit1_sessionindex + 1);

% remission_crit1_met = first_instance_under10(indices_less_than_10 ~= 0);