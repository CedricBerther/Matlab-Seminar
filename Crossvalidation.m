% Display prompt message in the command window
fprintf('Enter a number between 1 and 3 for the Dataset: ');

% Read input from the user
numberDataset = input('');

% Check if the input is a valid number between 1 and 3
while ~isnumeric(numberDataset) || numberDataset < 1 || numberDataset > 3
    fprintf('Invalid input. Please enter a number from 1 to 3: ');
    numberDataset = input('');
end

fileName = ['NewDemoData' num2str(numberDataset)];
load(fileName)

% Display the entered number
fprintf('You loaded the file: NewDemoData%d\n', numberDataset);

% Display prompt message in the command window
fprintf('How many folds do you want? ');

% Read input from the user
numberFolds = input('');

% Check if the input is a valid number higher than 1
while ~isnumeric(numberFolds) || numberFolds < 2 
    fprintf('Invalid input. Please enter a number higher than 1: ');
    numberFolds = input('');
end

% Display the entered number
fprintf('Your number of folds is: %d\n', numberFolds);


%% Scenario 2: Cross validation
% Example: Is their statistical evidence that we need a quadratic function
% to model the relationship between x and y?
h(1)=figure;

% Load data into x and y
x = data.x;
y = data.y;

% Make uniformly spaced x data, just for plotting
xPlot=linspace(min(x), max(x), 101);
N=numel(x);

% Plot x against y -- simple plot of the empirical data
subplot(1,2,1);
plot(x,y,'o')
xlabel('Intensity');
ylabel('Reaction time');

% Fit linear model to data:
[linearFitCoefficients]=polyfit(x,y,1);

% Evaluate the predictions of the linear model for the different data-points
yLinear=polyval(linearFitCoefficients,x);

% Calculate prediction accuracy by evaluating the norm of the difference between the actual y-values and the ones predicted 
% by the linear model:
mseLinear=norm(y-yLinear);

% Also evaluate linear model at plot-points for plotting:
yPlotLinear=polyval(linearFitCoefficients,xPlot);

% Repeat all steps above for the quadratic model:
[quadraticFitCoefficients]=polyfit(x,y,2);
yQuadratic=polyval(quadraticFitCoefficients,x);
mseQuadratic=norm(y-yQuadratic)
yPlotQuadratic=polyval(quadraticFitCoefficients,xPlot);

% Add the linear fit and quadratic fits to the empirical data:
hold on
plot(xPlot,yPlotLinear,'r','linewidth',2);
plot(xPlot,yPlotQuadratic,'g','linewidth',2);
legend('data','linear fit','quadratic fit');

%% Perform cross-validation:
% Check what type of cross-validation is run
if isnumeric(numberFolds)
    if numberFolds > N, error('Invalid crossvalidation parameter supplied: larger than number of datapoints!'), end
    if numberFolds < 2, error('Invalid crossvalidation parameter supplied: smaller than 2, the minimum number'), end
    disp(['Performing ' num2str(numberFolds) '-fold crossvalidation'])
elseif isstr(numberFolds)
    if strcmp(numberFolds, 'N')
        disp('Performing a leave-one-out cross-validation')
        numberFolds = N;
    else
        error('Invalid crossvalidation parameter supplied')
    end
end

% Calculate the step-size. Make sure it is an integer and warn the user  that the last fold will have a different size
% in case the number of datapoints cannot be divided by the number of folds without a remainder!
stepSize = floor(N/numberFolds);
remainder = rem(N, numberFolds);
disp(['N = ' num2str(N) ' datapoints and crossvalidation with ' ...
        num2str(numberFolds) '-folds:'])
if remainder > 0
    disp(['Note: There will be ' num2str(stepSize) ' datapoint(s) per fold in the test-set and ' num2str(N-stepSize) ' datapoint(s) in the training set' ...
       ' except in the last fold (' num2str(stepSize+remainder) ' test- and ' num2str(N-stepSize-remainder) ' training-datapoints)'])
else
    disp(['Note: There will be ' num2str(stepSize) ' datapoint(s) per fold in the test-set and ' num2str(N-stepSize) ' datapoint(s) in the training set'])
end

% Important: randomize the data once before starting the crossvalidation loop; shuffle both x and y data using the same
% randomisation!
randomIndex = randperm(length(x));
x = x(randomIndex);
y = y(randomIndex);

for k=1:numberFolds
   % Create the appropriate indices to select the correct number of points
   testIndex = [1+(k-1)*stepSize:k*stepSize];
   if k == numberFolds && remainder > 0
       testIndex = [testIndex testIndex(end)+[1:remainder]]; 
   end
   trainIndex = setdiff([1:N], testIndex);
   
   % Select the appropriate test and training sets in x and y
   xTrain = x(trainIndex);
   yTrain = y(trainIndex);
   xTest=x(testIndex);
   yTest=y(testIndex);
   
   % Fit model on training-set:
   [linearFitCoefficientsFold]=polyfit(xTrain,yTrain,1);
   [quadraticFitCoefficientsFold]=polyfit(xTrain,yTrain,2);

   % Calculate predictions on test-set:
   yLinearTest=polyval(linearFitCoefficientsFold,xTest);
   yQuadraticTest=polyval(quadraticFitCoefficientsFold,xTest);
    
   % Calculate error or quality of fit or goodness of fit on test-set:
   mseLinearFold(k)=norm(yTest-yLinearTest);
   mseQuadraticFold(k)=norm(yTest-yQuadraticTest);
end

subplot(1,2,2)
hist(mseLinearFold-mseQuadraticFold)
title('Linear minus quadratic MSE differences');

% Calculate the average MSE as well as percent of times the simpler model (linear) had a smaller test error than the
% more complex model (quadratic):
averageMSElinear=mean(mseLinearFold)
averageMSEquadratic=mean(mseQuadraticFold)
tmp = find(mseLinearFold <= mseQuadraticFold);
disp(['The linear model had a smaller test error on ' num2str(100*length(tmp)/numberFolds) '% of the ' ...
    'crossvalidation folds'])












