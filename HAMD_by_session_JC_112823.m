% % Transform Nili's Outcome and Static Predictor Table to Continuoius
% Process form % % 
clear all
%% load table
opts = delimitedTextImportOptions("NumVariables", 88);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["ID", "H_TOTALSCORE_0", "H_TOTALSCORE_1", "H_TOTALSCORE_2", "H_TOTALSCORE_3", "H_TOTALSCORE_4", "H_TOTALSCORE_5", "H_TOTALSCORE_6", "H_TOTALSCORE_7", "H_TOTALSCORE_8", "H_TOTALSCORE_9", "H_TOTALSCORE_10", "H_TOTALSCORE_11", "H_TOTALSCORE_12", "H_TOTALSCORE_13", "H_TOTALSCORE_14", "H_TOTALSCORE_15", "H_TOTALSCORE_16", "H_TOTALSCORE_17", "sg_session_n", "all_crit"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";

% Import the data
HAMD_by_session = readtable("/Users/jolinchou/Documents/PRIDE_Data_Analyses/id_HAMDtotal", opts);
HAMD_by_session_array = table2array(HAMD_by_session);
HAMD_by_session_array(:, 20:21) = [];

%% Detect sudden gains

% Initialize a logical array to store results
meetsCriteria = false(size(HAMD_by_session_array)); % Assuming 240 rows

% % Get the number of people in the table
% TotalSessions = 2:width(HAMD_by_session);

% Loop through each participant
for sessionindex = 3
    % Compute differences between consecutive scores for each participant
    HAMD_Delta =  HAMD_by_session_array(:, sessionindex + 1) - HAMD_by_session_array(:, sessionindex);

    % Calculate the percentage differences
    HAMD_percentdecrease =  (HAMD_by_session_array(:, sessionindex) - HAMD_by_session_array(:, sessionindex + 1)) ./ HAMD_by_session_array(:, sessionindex) * 100;

    % Extract prior 3 ECT sessions
    HAMD_pre = HAMD_by_session_array(:, (sessionindex - 2):(sessionindex));  

    if sessionindex == 3
         HAMD_pre = HAMD_by_session_array(:, (sessionindex - 1):(sessionindex));  
    end

    % Extract post 3 ECT sessions
    HAMD_post = HAMD_by_session_array(:, (sessionindex + 1):(sessionindex + 3));  

    % Calculate pre standard deviation 
    std_pre = nanstd((HAMD_pre), 0, 2);

    % Calculate post standard deviation 
    std_post = nanstd((HAMD_post), 0, 2);

    % Assuming HAMD_pre and HAMD_post are your arrays
    rows_with_NaN = any(isnan(HAMD_pre), 2) | any(isnan(HAMD_post), 2);

    % Assuming HAMD_pre and HAMD_post are your arrays
    rows_with_NaN_pre = sum(isnan(HAMD_pre), 2);
    rows_with_NaN_post = sum(isnan(HAMD_post), 2);

    % Rows with 0 NaN
    rows_with_zeros = (rows_with_NaN_pre == 0) & (rows_with_NaN_post == 0);

    % Rows with 1 NaN
    rows_with_1_NaN_pre = (rows_with_NaN_pre == 1) 

    | (rows_with_NaN_post == 1);

    % Rows with 1 NaN pre and post
    rows_with_1_NaN_prepost = (rows_with_NaN_pre == 1) & (rows_with_NaN_post == 1);

    % Rows with 2 NaNs pre or post
    rows_with_2_NaNs = (rows_with_NaN_pre == 2) | (rows_with_NaN_post == 2);

    % Indexing with rows_with_NaN to access rows with NaN values
    One_NaN_rows = find(rows_with_1_NaN);
    PrePost_NaN_rows = find(rows_with_1_NaN_prepost);
    bye_rows =  find(rows_with_2_NaNs);
    zero_rows = find(rows_with_zeros);

    doublecheck = length(One_NaN_rows) + length(PrePost_NaN_rows) + length(bye_rows) + length(zero_rows);

% % Iterate through rows
%     for rowIndex = 1:240
%         % Apply different conditions to specific rows
%         if any(isnan(HAMD_pre), 2) | any(isnan(HAMD_post), 2);
%             symptomfluctuation = 3.182 * sqrt(((HAMD_by_session_array(NaN_rows, sessionindex) - 1) .* std_pre.^2 + (HAMD_by_session_array(NaN_index, sessionindex + 1) - 1) .* std_post.^2) ./ (HAMD_by_session_array(NaN_index, sessionindex) + HAMD_by_session_array(NaN_index, sessionindex + 1) - 2));
%         else symptomfluctuation = 2.776 * sqrt(((HAMD_by_session_array(zero_index, sessionindex) - 1) .* std_pre.^2 + (HAMD_by_session_array(zero_index, sessionindex + 1) - 1) .* std_post.^2) ./ (HAMD_by_session_array(zero_index, sessionindex) + HAMD_by_session_array(zero_index, sessionindex + 1) - 2));
%         end
%     end
% 
    % Check if this row has NaN in either HAMD_pre or HAMD_post

        % for h = 1:length(3NaN_rows)
        %     3NaN_index = 3NaN_rows(j);
        %     symptomfluctuation =  * sqrt(((HAMD_by_session_array(2NaN_index, sessionindex) - 1) .* std_pre(2NaN_index).^2 + (HAMD_by_session_array(2NaN_index, sessionindex + 1) - 1) .* std_post(2NaN_index).^2) ./ (HAMD_by_session_array(2NaN_index, sessionindex) + HAMD_by_session_array(2NaN_index, sessionindex + 1) - 2));
        % 
        %     allsymptomfluctuation(NaN_index,:) = symptomfluctuation;
        % end 
        % 
        % for i = 1:length(2NaN_rows)
        %     2NaN_index = 2NaN_rows(j);
        %     symptomfluctuation = 4.303 * sqrt(((HAMD_by_session_array(2NaN_index, sessionindex) - 1) .* std_pre(2NaN_index).^2 + (HAMD_by_session_array(2NaN_index, sessionindex + 1) - 1) .* std_post(2NaN_index).^2) ./ (HAMD_by_session_array(2NaN_index, sessionindex) + HAMD_by_session_array(2NaN_index, sessionindex + 1) - 2));
        % 
        %     allsymptomfluctuation(NaN_index,:) = symptomfluctuation;
        % end 
        % 
        % for j = 1:length(1NaN_rows)
        %     1NaN_index = 1NaN_rows(j);
        %     symptomfluctuation = 3.182 * sqrt(((HAMD_by_session_array(1NaN_index, sessionindex) - 1) .* std_pre(1NaN_index).^2 + (HAMD_by_session_array(1NaN_index, sessionindex + 1) - 1) .* std_post(1NaN_index).^2) ./ (HAMD_by_session_array(1NaN_index, sessionindex) + HAMD_by_session_array(1NaN_index, sessionindex + 1) - 2));
        % 
        %     allsymptomfluctuation(NaN_index,:) = symptomfluctuation;
        % end 
        % 
        % for k = 1:length(zero_rows);
        %     zero_index = zero_rows(k)
        % symptomfluctuation = 2.776 * sqrt(((HAMD_by_session_array(zero_index, sessionindex) - 1) .* std_pre(NaN_index).^2 + (HAMD_by_session_array(zero_index, sessionindex + 1) - 1) .* std_post(NaN_index).^2) ./ (HAMD_by_session_array(zero_index, sessionindex) + HAMD_by_session_array(zero_index, sessionindex + 1) - 2));
        % 
        % allsymptomfluctuation(zero_index,:) = symptomfluctuation;
        % end

    % Criteria 3: Calculate Mpre - Mpost
    MeanDelta = mean(HAMD_pre, 2) - mean(HAMD_post, 2);

    % % Criteria 3: Calculate symptom fluctuation formula
    % symptomfluctuation = 2.776 * sqrt(((HAMD_by_session_array(:, sessionindex) - 1) .* std_pre.^2 + (HAMD_by_session_array(:, sessionindex + 1) - 1) .* std_post.^2) ./ (HAMD_by_session_array(:, sessionindex) + HAMD_by_session_array(:, sessionindex + 1) - 2));
    % 
    % Check if all three criteria are met for each score difference
    % criteriaCheck = HAMD_Delta >= -7 & HAMD_percentdecrease >= 25 &  MeanDelta > allsymptomfluctuation;

    meetsCriteria(:, sessionindex) = criteriaCheck;
end

% Identify rows that meet all criteria
rows_with_ones = any(meetsCriteria == 1, 2);
indices_with_ones = find(rows_with_ones);

 %% Initialize the number of rows and columns for 4x4 subplotx
numRows = 4;
numCols = 6;

% Calculate the total number of subplots
numSubplots = numRows * numCols;

% Create a cell array to hold the handles to individual figures
figureHandles = cell(1, numSubplots);

%% Generate figures for each participant

variableNames = {'sg_ppt'};

% Convert indicies_with_ones from a logical to a table format
rows_with_ones_table = array2table(rows_with_ones, 'variableNames', {'sg_ppt'});

% Concatenate sg_ppt to HAMD_by_session
HAMD_by_session = [HAMD_by_session, rows_with_ones_table];

% Get the number of people in the table
TotalParticipants = size(HAMD_by_session, 1);

% Loop through each person and create a figure
for personIndex = 1:TotalParticipants
    % Store all HAMD scores for the current person into an array
    individual_HAMD = HAMD_by_session(personIndex, 2:19);
    all_crit = HAMD_by_session.all_crit(personIndex);
    sg_session_n = HAMD_by_session.sg_session_n(personIndex);
    sg_ppt = HAMD_by_session.sg_ppt(personIndex);

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
    if sg_ppt == 1
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-r*')

        % Add labels, titles, legends, etc. to the figure
        title(['Subject ', num2str(ID), ' SG: session ', num2str(sg_session_n)]);
        xlabel('ECT Treatment Visit');
        ylabel('HAM-D Score');

        %Set xlim and ylim
        xlim([0, max(18)]);
        ylim([min(0), max(60)]);
    else

    plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-*')

    % Add labels, titles, legends, etc. to the figure
    title(['Subject ', num2str(ID)]);
    xlabel('ECT Treatment Visit');
    ylabel('HAM-D Score');

    %Set xlim and ylim
    xlim([0, max(18)]);
    ylim([min(0), max(60)]);

    end
end


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
