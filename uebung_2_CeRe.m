%%  1. Clear all variables in the workspace

clear;

%%  2. Find out which directory you are currently in, and display the names 
%      of all files in the current directory

disp(pwd);
disp(ls);

%%  3. Without opening any of the three files, find out which of them contains a 
%      variable named BigRandomMatrix. Load only this variable from the file

DirList = dir(fullfile(pwd, '*.mat'));
Data = cell(1, length(DirList));
for k = 1:length(DirList)
  test = ismember('BigRandomMatrix', who('-file', DirList(k).name));
  if(test == 1)
      disp("Variable ist in Datei: " + DirList(k).name);
      variable = load(DirList(k).name,"BigRandomMatrix");
  end
end

%%  4. Clear all variables; load TestFile1.mat, and display the variables in the
%      workspace. Which variable takes the most amount of memory? 
%      Clear this variable. Save the remaining variables under TestFile1Small.mat
clear;
load('TestFile1.mat');
whos;
clear d;
save('TestFile1Small.mat')



%%  5. Save the variables a b c d from TestFile1.mat in a file called 
%      TestFile1Letters.mat

clear;
load('TestFile1.mat');
clear BigRandomMatrix;
save('TestFile1Letters.mat')

%%  6. Type in numel someText and numel(someText). Why do the two function
%      calls lead to different results?


%%  7. Display the contents of the variable someText. Load the file TestFile2.mat,
%      and display the contents of someText again, ditto for TestFile3.mat

%%  8. Now, append the variables someText and logicalEye to TestFile1Small.mat.
%      Did MATLAB overwrite the variable someText (which already exists in
%      TestFile1Small.mat) or not?
%%  9. Clear the workspace, and load TestFile1.mat and TestFile2.mat into
%      structs. Now, you can display the content of the two someText variables
%      simultaneously
%%  10. Write an m-file which takes as input a filename, checks whether the file
%       exists, loads the file if yes, and displays a warning otherwise


