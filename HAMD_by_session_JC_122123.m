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
 allHAMD_Delta = zeros(size(HAMD_by_session_array)); 
 sg_reversalvalue = [];

% % Get the number of people in the table
% TotalSessions = 2:width(HAMD_by_session);

% Loop through each participant
for sessionindex = 3:17;
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
    MeanDelta = nanmean(HAMD_pre, 2) - nanmean(HAMD_post, 2);
    
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

variableNames = {'sg_ppt'};

% Convert indicies_with_ones from a logical to a table format
rows_with_ones_table = array2table(rows_with_ones_sg_reversed, 'variableNames', {'sg_ppt'});

% Concatenate sg_ppt to HAMD_by_session
HAMD_by_session = [HAMD_by_session, rows_with_ones_table];

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


fprintf('Our total sudden gainers: %d\n', our_sgsum);
fprintf('Nili total sudden gainers: %d\n', Nili_sgsum);
fprintf('Total matching sudden gainers: %d\n', sum_matching_gainers);

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
    if sg_ppt == 1 & all_crit == 1
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-m*')

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

        
    elseif sg_ppt == 1 & all_crit == 0 
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
        
    elseif sg_ppt == 0 & all_crit == 1
        plot(0:(width(individual_HAMD) - 1), table2array(individual_HAMD), '-b*')

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
