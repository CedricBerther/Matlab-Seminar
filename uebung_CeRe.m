function [ out ] = uebung_CeRe( myFileName, myDataDirectory)
%  AnalyseData_CeRe  reads your matlab-file and gives you a 
%  theoretical and your empirical data A and B, that you loaded from disk as a plot.
%  Y = AnalyseData_CeRe(myFileName, myDataDirectory) will use the
%  path that you send as a string.
%
%  If you dont want to write the name and path of your file by hand, you can
%  also send no variables with the function, this way you can choose your
%  file via UI. 
%% Define variables
TITLE_XLABLE = 'x-values';
TITLE_YLABLE = 'y-values';
DISP_PATH_NOT_EXISTING = 'Path does not exist';
DISP_IS_FILE = 'The specified file does not exist in the directory';
DISP_EMPTY = 'Variables are not allowed to be empty';
TITLE_UIGETFILE = 'Select a MATLAB .mat data-file';
DISP_LOADED_FILE = 'Loaded the following file: ';
DISP_CANCEL = 'User selected Cancel';
SYMBOLLIST = ['o', 'x', 's', 'd','.', '^', 'v', '>', '<', '*','p','h','+'];
SYMBOL = @(k) sprintf('r%s-', SYMBOLLIST(mod(k-1,length(SYMBOLLIST))+1));
COLORSTRING = 'kbmryg';
TEXT_A_MEAN = 'Difference of A in mean = ';
TEXT_B_MEAN = 'Difference of B in mean = ';
TEXT_A_STD = 'Difference of A in std = ';
TEXT_B_STD = 'Difference of B in std = ';
%% See how many variables where given
%  Use file when filename and path is given and 
%  see if Path exists, file is real and if given variables are not empty.
%  Use UI when nothing or filename is given.
switch nargin
    case 2
        assert(isfolder(myDataDirectory), DISP_PATH_NOT_EXISTING)
        assert(isfile(fullfile(myDataDirectory, myFileName)), DISP_IS_FILE)
        assert(~isempty(myFileName) && ~isempty(myDataDirectory), DISP_EMPTY)
    otherwise
        [myFileName, myDataDirectory] = uigetfile('*.mat', TITLE_UIGETFILE);
        if isequal(myFileName,0)
            error(DISP_CANCEL);
        else
            disp([DISP_LOADED_FILE fullfile(myDataDirectory, myFileName)]);
        end
end
%%  Load data and fill variables
%   Load data from file
completeAccessPath = fullfile(myDataDirectory, myFileName);
data = load(completeAccessPath);
%   Load empirical data from file
A = data.A;
B = data.B;
x = data.x;
%   Calculate theoretical functions 
yA = 3 * exp(-x) - 1;
yB = x.^2 - 1;
%   Calculate mean & standard deviation difference  
meanDifferenceA = norm(mean(yA - A));
meanDifferenceB = norm(mean(yB - B));
stdDifferenceA = norm(std(yA - A));
stdDifferenceB = norm(std(yB - B));
%  Smooth out theoretical lines. 
xSmooth = linspace(min(x), max(x), 200);
yASmooth = spline(x,yA,xSmooth);
yBSmooth = spline(x,yB,xSmooth);
%   Create text for output 
textAmean = [TEXT_A_MEAN num2str(meanDifferenceA)];
textBmean = [TEXT_B_MEAN num2str(meanDifferenceB)];
textAstd = [TEXT_A_STD num2str(stdDifferenceA)];
textBstd = [TEXT_B_STD num2str(stdDifferenceB)];
%%  Plot Data
%   Use different symbols and color for each plot, plot A, B and the
%   theoretical functions. Write text with mean and std in plot.
figure('name',myFileName,'NumberTitle','off');
plot(x,A,SYMBOL(1), 'Color', COLORSTRING(1));
hold on 
plot(x,B,SYMBOL(2), 'Color', COLORSTRING(2));
hold on 
plot(xSmooth,yASmooth,'Color', COLORSTRING(3));
hold on 
plot(xSmooth,yBSmooth, 'Color', COLORSTRING(4));
hold off
text(0.1,2.8,{textAmean, textBmean, textAstd, textBstd},'Color','red','FontSize',10)
title(['Plot of ' myFileName]);
xlabel(TITLE_XLABLE)
ylabel(TITLE_YLABLE)
%%  Display mean and std in command window
disp([textAmean newline textBmean newline textAstd newline textBstd]);
%%  If user wants output -> give back mean & std
switch nargout
    case 1 
        out = [meanDifferenceA stdDifferenceA meanDifferenceB stdDifferenceB];
end
end
