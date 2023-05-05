function [out] = CrossvalidationFunction(numberOfFolds, dataSet, B)
% [out] = CrossvalidationFunction(numberOfFolds, dataSet, B) 
%
% Modified flexibleCrossvalidationScript to be a function.
% Three input parameters:
%   1. numberOfFolds = what type of crossvalidation to perform
%   2. dataSet = number 1, 2 or 3 to load NewDemoData1, 2 or 3
%   3. B = number of repetitions of crossvalidation to perform
%       (Note: if numberOfFolds is equal to "N", i.e. leave-one-out
%       crossvalidation, then B is set to 1, as repeating this is not
%       useful.)
%
% FAW       26.04.2018  First version
%           30.04.2023  Clean-up and more feedback in Command Window; 


%% Set defaults, load file, get data, initialize data struct
VERBOSE = 'true';                        % displays results in the console and shows figure
NUMBER_X_SMOOTH = 200;

if nargin < 1, numberOfFolds = 10; end
if nargin < 2, dataSet = 2; end
if nargin < 3, B = 2; end
fileName = ['NewDemoData' num2str(dataSet)];
if strcmp(VERBOSE, 'true')
    disp(' ')
    disp(['Loading file named "' fileName '" from directory ' cd])
end
load(fileName);

x = data.x;
y = data.y;
N=numel(x);
xPlot=linspace(min(x), max(x), NUMBER_X_SMOOTH);

out.empiricalData.x = x;
out.empiricalData.y = y;
out.empiricalData.N = N;
out.empiricalData.fileName = fileName;
out.empiricalData.directory = cd;


%% Fit a linear and quadratic function to the data
%Fit linear model to data:
[linearFitCoefficients]=polyfit(x,y,1);
%evaluate the predictions of the lineaer model for the different
%data-points
yLinear=polyval(linearFitCoefficients,x);
%calculate prediction accuracy by evaluating the norm of the difference
%between the actual y-values and the ones predicted by the linear model:
mseLinear=norm(y-yLinear);
%now, also evaluate linear model at plot-points for plotting:
yPlotLinear=polyval(linearFitCoefficients,xPlot);

%repeat the same steps for the quadratic model:
[quadraticFitCoefficients]=polyfit(x,y,2);
yQuadratic=polyval(quadraticFitCoefficients,x);
mseQuadratic=norm(y-yQuadratic);
yPlotQuadratic=polyval(quadraticFitCoefficients,xPlot);

if strcmp(VERBOSE, 'true')
   disp(' ')
   disp('Whenever you fit a more complex model (typically with more free parameters) the error of the fit is lower')
   disp(['MSE of the linear model:    ' num2str(mseLinear)])
   disp(['MSE of the quadratic model: ' num2str(mseQuadratic)])
end
out.empiricalData.linearMSE = mseLinear;
out.empiricalData.quadraticMSE = mseQuadratic;
out.empiricalData.linearParams = linearFitCoefficients;
out.empiricalData.quadraticParams = quadraticFitCoefficients;

%% Plot raw data and best fitting linear and quadratic functions
figure(1), clf, hold off
subplot(1,2,1);
plot(x,y,'o')
xlabel('Intensity');
ylabel('Reaction time');
% plot linear fit to the data and quadratic fit into our plot:
hold on
plot(xPlot,yPlotLinear,'r','linewidth',2);
plot(xPlot,yPlotQuadratic,'g','linewidth',2);
legend('data','linear fit','quadratic fit');

%% Prepare cross-validation
% check what type of cross-validation is run
if isnumeric(numberOfFolds)
    if numberOfFolds > N, error('Invalid crossvalidation parameter supplied: larger than number of datapoints!'), end
    if numberOfFolds < 2, error('Invalid crossvalidation parameter supplied: smaller than 2, the minimum number'), end
    disp(['Performing ' num2str(numberOfFolds) '-fold crossvalidation'])
elseif isstr(numberOfFolds)
    if strcmp(numberOfFolds, 'N')
        disp('Performing a leave-one-out cross-validation')
        numberOfFolds = N;
        if B ~= 1               % If leave-one-out, set B to 1 if not at 1
           B = 1;
           disp('No need to repeat a leave-one-out cross-validation; thus parameter B set to 1')
        end
    else
        error('Invalid crossvalidation parameter supplied')
    end
end
out.cv.numberOfFolds = numberOfFolds;
out.cv.B = B;

% calculate the step-size. Make sure it is an integer and warn the user
% that the last fold will have a different size in case the number of
% datapoints cannot be divided by the number of folds without a remainder!
stepSize = floor(N/numberOfFolds);
out.cv.testSetSize = stepSize;
out.cv.trainingSetSize = N - stepSize;

remainder = rem(N, numberOfFolds);
if strcmp(VERBOSE, 'true')
    disp(' ')
    disp(['N = ' num2str(N) ' datapoints and crossvalidation with ' ...
            num2str(numberOfFolds) '-folds:'])
    if remainder > 0
        disp(['Note: There will be ' num2str(stepSize) ' datapoint(s) per fold in the test-set and ' num2str(N-stepSize) ' datapoint(s) in the training set' ...
           ' except in the last fold (' num2str(stepSize+remainder) ' test- and ' num2str(N-stepSize-remainder) ' training-datapoints)'])
        out.cv.testSetSizeLastFold = stepSize+remainder;
        out.cv.trainingSetSizeLastFold = N - stepSize-remainder;

    else
        disp(['Note: There will be ' num2str(stepSize) ' datapoint(s) per fold in the test-set and ' num2str(N-stepSize) ' datapoint(s) in the training set'])
    end
end


%% Perform B-times a cross-validation:
out.cv.linearParams = zeros(B*numberOfFolds, 2);
out.cv.quadraticParams = zeros(B*numberOfFolds, 3);
out.cv.linearMSE = zeros(B, numberOfFolds);
out.cv.quadraticMSE = zeros(B, numberOfFolds);
out.cv.averageLinearMSE = zeros(B,1);
out.cv.averageQuadraticMSE = zeros(B,1);
for j = 1:B
    % important: randomize the data once before starting the crossvalidation
    % loop on very of the B iterations;
    % Shuffle both x and y data using the same randomisation as the data
    % comes in pairs!
    randomIndex = randperm(length(x));
    x = x(randomIndex);
    y = y(randomIndex);
    
    for k=1:numberOfFolds
       % create the appropriate indices to select the correct number of
       % points.
       testIndex = [1+(k-1)*stepSize:k*stepSize];
       % In case there is a remainder, i.e. the test- and training set
       % sizes in the last fold are different, we have to add the remainder
       % indices to the testIndex:
       if k == numberOfFolds && remainder > 0
           testIndex = [testIndex testIndex(end)+[1:remainder]]; 
       end
       trainIndex = setdiff([1:N], testIndex);

       % select the appropriate test and training sets in x and y
       xTrain = x(trainIndex);
       yTrain = y(trainIndex);
       xTest=x(testIndex);
       yTest=y(testIndex);

       %fit model on training-set:
       [linearFitCoefficientsFold]=polyfit(xTrain,yTrain,1);
       [quadraticFitCoefficientsFold]=polyfit(xTrain,yTrain,2);

       %calculate predictions on test-set:
       yLinearTest=polyval(linearFitCoefficientsFold,xTest);
       yQuadraticTest=polyval(quadraticFitCoefficientsFold,xTest);

       %calculate goodness of fit on test-set
       out.cv.linearMSE(j, k)=norm(yTest-yLinearTest);
       out.cv.quadraticMSE(j, k)=norm(yTest-yQuadraticTest);
    end

    %average mse:
    out.cv.averageLinearMSE(j) = mean(out.cv.linearMSE(j, :) );
    out.cv.averageQuadraticMSE(j) = mean(out.cv.quadraticMSE(j, :) );
end

% Plot the histogram of the error
subplot(1,2,2)
linearMSEs = reshape(out.cv.linearMSE, [], 1);
qudraticMSEs = reshape(out.cv.quadraticMSE, [], 1);
hist(linearMSEs - qudraticMSEs)
title({'Differences in MSE on test sets';'Positive difference: MSE_{x} > MSE_{x^2}'});


fitError = linearMSEs - qudraticMSEs;
tmp = find(fitError > 0);
out.cv.percentageTestErrorLinearLarger = 100*length(tmp)/length(fitError);


% More feedback in the Command Window
disp(' ')
disp(['The linear model had a smaller test error on ' num2str(100-out.cv.percentageTestErrorLinearLarger) '% of the ' ...
    'crossvalidation folds, given B = ' num2str(B) ' repeats.'])
disp('Average test errors across all folds and repeats were:')
disp(['Linear    model: ' num2str(mean(out.cv.averageLinearMSE))])
disp(['Quadratic model: ' num2str(mean(out.cv.averageQuadraticMSE))])


end









