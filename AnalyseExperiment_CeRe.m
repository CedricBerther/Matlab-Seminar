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
% Create a plot with different colors for each session
figure;

% Plot Figure 1
subplot(1, 2, 1);  % Crete a subplot grid of 1 row and 2 columns, select the first subplot
hold on;
scatter(Level_1 - 0.15, RT_1, COLOR(1), SYMBOL(1));
scatter(Level_2 + 0.15, RT_2, COLOR(2), SYMBOL(2));
hold off;

title('Figure 1');
xlabel(FIGURE_1_X_LABEL);
ylabel(FIGURE_1_Y_LABEL);
legend(date_1, date_2, 'Location', 'north');

%% Plot Figure 2
%HIER WEITERMACHEN. RMSE FÜR ZWEI LEVEL WURDEN BERECHNET
xPlot = linspace(min(Level_2), max(Level_2), NUMBER_X_SMOOTH);

[FitCoefficients_1] = polyfit(Level_2, RT_2, 1);
y_1 = polyval(FitCoefficients_1, Level_2);
yPlot_1 = polyval(FitCoefficients_1, xPlot);
A = rmse(y_1, RT_2);
disp(A);

[FitCoefficients_2]=polyfit(Level_2,RT_2,2);
yPlot_2=polyval(FitCoefficients_2,xPlot);

[FitCoefficients_3]=polyfit(Level_2,RT_2,3);
yPlot_3=polyval(FitCoefficients_3,xPlot);

[FitCoefficients_4]=polyfit(Level_2,RT_2,4);
yPlot_4=polyval(FitCoefficients_4,xPlot);

[FitCoefficients_5]=polyfit(Level_2,RT_2,5);
y_5 = polyval(FitCoefficients_5, Level_2);
yPlot_5=polyval(FitCoefficients_5,xPlot);
B = rmse(y_5, RT_2);
disp(B);

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
legend(date_2, 'Location', 'north');

%% Add parameters to struct

%% Create File

%% Leave-one-out cross-validation

%% Add information to struct

%% Plot Figure 3

%% Non-parametric bootstrap test 

%% Add information about benefit of training to struct

%% Save data

%% If user wants output -> give back output
switch nargout
    case 1 
        out = ["OUTPUT"];
end
end