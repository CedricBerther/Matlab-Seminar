function [] = uebung_CeRe()
%%  Load file into workplace

teststruct = load('teststruct.mat');


%%  Remove the fields xeps and yrate from the struct seq

teststruct.seq = rmfield(teststruct.seq,"xeps");
teststruct.seq = rmfield(teststruct.seq,"yrate");

%%  Make a vector which contains all the entries of seq.T

 vector_t = [teststruct.seq.T];
 
%%  Make a struct (using structfun) which extracts the maximum element from each field in seq(1)

max_seq1 = structfun(@max,teststruct.seq(1), UniformOutput=false);

%%  Extract (using structfun) the size of each entry of seq(1)

size_seq1 = structfun(@size,teststruct.seq(1), UniformOutput=false);

%%  Make a cell which contains each entry of seq.x

cell_x = {teststruct.seq.x};

%%  Using cellfun, find out the number of elements of each of the members of the cell-array you generated in the previous exercise

n_of_elements_Cell = cellfun(@numel, cell_x);

%%  Find out for which indices seq.T does NOT match the number of elements of seq.x

find_indices = find(vector_t ~= n_of_elements_Cell);
find_indices

%%  Convert seq to a cell-array using struct2cell

cell_seq = struct2cell(teststruct.seq);

end
