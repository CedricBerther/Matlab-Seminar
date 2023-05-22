function [ out ] = AnalyseExperiment_CeRe()
%  Explanation

%% Define variables
FILENAME = "ProblemSolvingExperiment.mat";
UI_TITLE = 'Select a MATLAB .mat data-file';
UI_DISP_LOADED_FILE = 'Loaded the following file: ';
UI_DISP_CANCEL = 'User selected Cancel';
UI_DISP_COMPLETED = 'User selection completed';
DATADIRECTORY = pwd;
FIGURE_1_TITLE = "Figure 1";
FIGURE_1_X_LABEL = "Level";
FIGURE_1_Y_LABEL = "Reactiontime";
SYMBOLLIST = ['o', 'x', 's', 'd','.', '^', 'v', '>', '<', '*','p','h','+'];
SYMBOL = @(k) sprintf('%s', SYMBOLLIST(mod(k-1,length(SYMBOLLIST))+1));
COLORLIST = ['r', 'b', 'g', 'm', 'c'];
COLOR = @(k) sprintf('%s', COLORLIST(mod(k-1,length(COLORLIST))+1));
NUMBER_X_SMOOTH = 200;
NUMBER_FOLDS = 10;
%% Find File

if exist(FILENAME, 'file') ~= 2
    [FILENAME, DATADIRECTORY] = uigetfile('*.mat', UI_TITLE);
    if isequal(FILENAME,0)
        error(UI_DISP_CANCEL);
    else
        disp(UI_DISP_COMPLETED);
    end
end
completeAccessPath = fullfile(DATADIRECTORY, FILENAME);
data = load(completeAccessPath);
disp([UI_DISP_LOADED_FILE fullfile(completeAccessPath)]);

%% Load data and fill variables

structA = data.a;
structA.path = completeAccessPath;

data_1 = structA.session(1).data; 
data_2 = structA.session(2).data;
date_1 = datestr(structA.session(1).date); 
date_2 = datestr(structA.session(2).date);

RT_1 = data_1(:, 1);
Level_1 = data_1(:, 2);
RT_2 = data_2(:, 1);
Level_2 = data_2(:, 2);

%% Plot Figure 1

figure;
subplot(1, 2, 1); 
hold on;
scatter(Level_1 - 0.15, RT_1, COLOR(1), SYMBOL(1));
scatter(Level_2 + 0.15, RT_2, COLOR(2), SYMBOL(2));
hold off;

title('Figure 1');
xlabel(FIGURE_1_X_LABEL);
ylabel(FIGURE_1_Y_LABEL);
legend(date_1, date_2, 'Location', 'north');

%% Plot Figure 2

xPlot = linspace(min(Level_2), max(Level_2), NUMBER_X_SMOOTH);

[FitCoefficients_1] = polyfit(Level_2, RT_2, 1);
y_1 = polyval(FitCoefficients_1, Level_2);
yPlot_1 = polyval(FitCoefficients_1, xPlot);
rmse_1 = rmse(y_1, RT_2);

[FitCoefficients_2]=polyfit(Level_2,RT_2,2);
y_2 = polyval(FitCoefficients_2, Level_2);
yPlot_2=polyval(FitCoefficients_2,xPlot);
rmse_2 = rmse(y_2, RT_2);

[FitCoefficients_3]=polyfit(Level_2,RT_2,3);
y_3 = polyval(FitCoefficients_3, Level_2);
yPlot_3=polyval(FitCoefficients_3,xPlot);
rmse_3 = rmse(y_3, RT_2);

[FitCoefficients_4]=polyfit(Level_2,RT_2,4);
y_4 = polyval(FitCoefficients_4, Level_2);
yPlot_4=polyval(FitCoefficients_4,xPlot);
rmse_4 = rmse(y_4, RT_2);

[FitCoefficients_5]=polyfit(Level_2,RT_2,5);
y_5 = polyval(FitCoefficients_5, Level_2);
yPlot_5=polyval(FitCoefficients_5,xPlot);
rmse_5 = rmse(y_5, RT_2);

subplot(1, 2, 2);  % Select the second subplot
hold on;
scatter(Level_2 + 0.15, RT_2, COLOR(2), SYMBOL(2));
plot(xPlot,yPlot_1, COLOR(1),'linewidth',1);
plot(xPlot,yPlot_2, COLOR(2),'linewidth',1);
plot(xPlot,yPlot_3, COLOR(3),'linewidth',1);
plot(xPlot,yPlot_4, COLOR(4),'linewidth',1);
plot(xPlot,yPlot_5, COLOR(5),'linewidth',1);
hold off;

title('Figure 2');
xlabel(FIGURE_1_X_LABEL);
ylabel(FIGURE_1_Y_LABEL);
text(0.8, 24, sprintf('RMSE 1: %.2f', rmse_1), 'FontSize', 10, 'FontWeight', 'bold');
text(0.8, 22, sprintf('RMSE 2: %.2f', rmse_2), 'FontSize', 10, 'FontWeight', 'bold');
text(0.8, 20, sprintf('RMSE 3: %.2f', rmse_3), 'FontSize', 10, 'FontWeight', 'bold');
text(0.8, 18, sprintf('RMSE 4: %.2f', rmse_4), 'FontSize', 10, 'FontWeight', 'bold');
text(0.8, 16, sprintf('RMSE 5: %.2f', rmse_5), 'FontSize', 10, 'FontWeight', 'bold');

%% Add parameters to struct

structA.parameter_fit_1 = y_1;
structA.parameter_fit_2 = y_2;
structA.parameter_fit_3 = y_3;
structA.parameter_fit_4 = y_4;
structA.parameter_fit_5 = y_5;

structA.rmse_1 = rmse_1;
structA.rmse_2 = rmse_2;
structA.rmse_3 = rmse_3;
structA.rmse_4 = rmse_4;
structA.rmse_5 = rmse_5;

%% Save File

data.a = structA;
save('ProblemSolvingExperiment_2.mat', '-struct', 'data');

%% Leave-one-out cross-validation

% Important: randomize the data once before starting the crossvalidation loop; shuffle both x and y data using the same
% randomisation!

N = numel(RT_2);
stepSize = floor(N/NUMBER_FOLDS);
remainder = rem(N, NUMBER_FOLDS);

rmsea = zeros(1, 5);

for j = 1:5
    % important: randomize the data once before starting the crossvalidation
    % loop on very of the B iterations;
    % Shuffle both x and y data using the same randomisation as the data
    % comes in pairs!
    randomIndex = randperm(length(RT_2));
    x_2 = RT_2(randomIndex);
    y_2 = Level_2(randomIndex);
    
    for k = 1:NUMBER_FOLDS
        % create the appropriate indices to select the correct number of
        % points.
        testIndex = [1 + (k-1) * stepSize:k * stepSize];
        % In case there is a remainder, i.e. the test- and training set
        % sizes in the last fold are different, we have to add the remainder
        % indices to the testIndex:
        if k == NUMBER_FOLDS && remainder > 0
            testIndex = [testIndex testIndex(end) + [1:remainder]]; 
        end
        trainIndex = setdiff([1:N], testIndex);

        % select the appropriate test and training sets in x and y
        xTrain = x_2(trainIndex);
        yTrain = y_2(trainIndex);
        xTest = x_2(testIndex);
        yTest = y_2(testIndex);

        % fit model on training-set:
        [linearFitCoefficientsFold] = polyfit(xTrain, yTrain, j);
        yLinearTest = polyval(linearFitCoefficientsFold, xTest);

        % calculate RMSEA on test-set
        linearRMSEA(j, k) = sqrt(mean((yTest - yLinearTest).^2));
    end

    % average RMSEA:
    averageLinearRMSEA(j) = mean(linearRMSEA(j, :));

    % store the RMSEA values for each polynomial degree
    rmsea(j) = averageLinearRMSEA(j);
end

% find the polynomial degree with the lowest average RMSEA
[~, bestPolynomial] = min(rmsea);


%% Add information to struct
structA.bestPolynomial = bestPolynomial;
%% Plot Figure 3

%% Non-parametric bootstrap test 

    session1_rt_difficulty_1 = data_1(data_1(:, 2) == 1, 1);
    session2_rt_difficulty_1 = data_2(data_2(:, 2) == 1, 1);

    session1_rt_difficulty_10 = data_1(data_1(:, 2) == 10, 1);
    session2_rt_difficulty_10 = data_2(data_2(:, 2) == 10, 1);

    N=numel(session1_rt_difficulty_10);
    numBootstraps = 5000; 
    bootstrapMeans = zeros(numBootstraps, 1);
    
    for i = 1:numBootstraps
        bootstrapIndex=randi(N,N,1);

        xbootstrapped=session1_rt_difficulty_10(bootstrapIndex);
        ybootstrapped=session2_rt_difficulty_10(bootstrapIndex);
    
        bootstrapMeans(i) = mean(ybootstrapped) - mean(xbootstrapped);

    end

    confidence_level = 0.95;  % Confidence level
    lower_percentile = (1 - confidence_level) / 2;
    upper_percentile = 1 - lower_percentile;

    bootstrap_mean_diff_sorted = sort(bootstrapMeans);
    lower_bound = bootstrap_mean_diff_sorted(round(lower_percentile * numBootstraps));
    upper_bound = bootstrap_mean_diff_sorted(round(upper_percentile * numBootstraps));

    if lower_bound <= 0 && upper_bound >= 0
        disp('The 95% confidence interval contains zero for task difficulty 1.');
    else
        disp('The 95% confidence interval does not contain zero for task difficulty 1.');
    end
    
    if lower_bound <= 0 && upper_bound >= 0
        disp('The 95% confidence interval contains zero for task difficulty 10.');
    else
        disp('The 95% confidence interval does not contain zero for task difficulty 10.');
    end

%% Add information about benefit of training to struct

p_values = zeros(10, 1);

for difficulty = 1:10
    session1_rt_difficulty = data_1(data_1(:, 2) == difficulty, 1);
    session2_rt_difficulty = data_2(data_2(:, 2) == difficulty, 1);

    [h,p] = ttest2(session1_rt_difficulty,session2_rt_difficulty)
    
    p_values(difficulty) = p;
end

significant_difficulty = find(p_values < 0.05, 1);

if isempty(significant_difficulty)
    disp('No difficulty level has a statistically significant difference.');
else
    disp(['The lowest difficulty level with a statistically significant difference is ', num2str(significant_difficulty)]);
end

structA.lowest_level_significant = significant_difficulty;

%% Save data
data.a = structA;
save('ProblemSolvingExperiment_2.mat', '-struct', 'data');
%% If user wants output -> give back output
switch nargout
    case 1 
        out = ["OUTPUT"];
end
end