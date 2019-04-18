function [ outMat,minMatIdx ] = hammDist(truncatedPep, fileName, hammdist)
% % Written by Leah V. Sibener 20141001
% % Edited by Timothy F. Miles 20190114

% This function will output a nx3 matrix of paired peptides and their reverse
% hamming distance.  In this version of the program the hamming distance is
% computed and then subtracted from the length of the peptide.  This will
% make it easier to use cytoscape to make the width of the lines mean
% more related.

% Input for this function is a cell array of peptides (you will import with
% directions below), a file name you wish to export it to (e.g.
% 'Test_hammingDistances'), and the minimum acceptable number of identical
% positions between pairs of sequences. 

%   Detailed Explanation: Align your sequences, take the variable region of interest 
%   and copy and paste it into an excel file.
      
%   Make sure that this excel file is in the same folder as this program.
%   OR ELSE IT WILL NOT WORK.

%   When you use the import wizard (just double click on excel file on
%   left) import all of the cells (peptides in column) as a CELL ARRAY.
%   Under text options, make sure 'CELL ARRAY with character vectors' is
%   selected.

%   Note: all of these peptides must be the same length or
%   else you cannot compute the hamming distance.

% HERE IS AN EXAMPLE OF HOW TO EXECUTE THE FUNCTION:
%   after you have imported a cell array lets say you call it Test1
%   hammingDist(Test1, 'Test1_hamming', 5) this will then be
%   exported to a txt file called Test1_hamming with five or greater 
%   reverse hamming distance.

N = length(truncatedPep);

% pre-allocation of where hamming matrix will go
hammingMat = zeros(N,N);

% THIS LOOP gets hamming distance
% length of peptide is M
M = length(truncatedPep{1,1});
for i= 1:N-1 % all peptides
    for j = i+1:N % compares to peptides but no duplicates (makes diagonal matrix)
        str1=truncatedPep{i}; %eventually all the peptides
        str2=truncatedPep{j}; %compare to peptide j
        sum1=0; %before comparing they have nothing the same
        for k=1:length(str1) % for the length of the peptide
            sum1=sum1+(str1(k)==str2(k)); % evaluates number of same and adds it to zero
        end
        hammingMat(i,j) = length(str1)-sum1; %this is a diagonal matrix
    end
end

% Make reflection to make the whole matrix
hammingMatreflect = hammingMat';

% full matrix is an N by N matrix of all the pairwise hamming distances
fullMatrix = hammingMat+hammingMatreflect; %addition of matricies
fullMatrix = fullMatrix+(eye(N,N)*10^6); % along the diagonal i have made it equal to 10^6, otherwise zeros along the diagonal
% would make finding the minimum hard

%find minimum values for each peptide by row
for i = 1:N
    [value] = min(fullMatrix(i,:),[],2);
    minMatValues(i) = value;
    minMatIdx(i,:) = value==fullMatrix(i,:);
   
end

% make list for output file
k = 1; % counter
lengthPep = length(truncatedPep{i});
 
for i = 1:N
    idxRow = minMatIdx(i,:);
    listIdx = find(idxRow);
    for j = listIdx
        if lengthPep-fullMatrix(i,j)>= hammdist
        outMat{k,1} = truncatedPep{i}; %all the peptides in first column
        outMat{k,2} = truncatedPep{j}; %minimum hamming distance in the second
        outMat{k,3} = lengthPep-fullMatrix(i,j); %length-hamming dist. in third
        k = k+1;
    end
  end


hamPairs = cell(outMat(cell2mat(outMat(:,3)) >= hammdist, :)); %data enter as cell array
pairLength = length(hamPairs); %length of array

for i = 1:pairLength
    seq2Trim1 = hamPairs{i,1};
    hamPairs{i,1} = seq2Trim1(1:lengthPep);
    
    seq2Trim2 = hamPairs{i,2}; %second column
    hamPairs{i,2} = seq2Trim2(1:lengthPep);
    
end

redVect = zeros(pairLength, 1); %stores vector of where redundant sequences are

count = 0; %counts non redundent
countRedundent = 0; %counts redundent


for i = 1:pairLength
  
   seq1 = hamPairs{i,1}; %column 1 seq
   seq2 = hamPairs{i,2}; %column 2 seq

   
   for j = 1:pairLength  % compare to all other pairs
       if j ~= i % unless its the same row
       seq3 = hamPairs{j,2}; % column 2 seq
       seq4 = hamPairs{j,1}; % column 1 seq
       
       tf1 = strcmp(seq1, seq3); %compare seq1 to seq3 etc etc
       tf2 = strcmp(seq2, seq4);
       tf3 = strcmp(seq1, seq4);
       tf4 = strcmp(seq2, seq3);
       
       if (tf1 && tf2) || (tf3 && tf4) == 1 % if intereaction is redundent

           redVect(i) = 1;  % place 1 evertime there is a redundant pair of sequences
           hamPairs{j,1} = ''; % replace second instance with no string
          hamPairs{j,2} = ''; 
          hamPairs{j,3} = 0; %replace hamming dist with 0
       end
       end
   end
end

% now collating all nonredundnet hamming distance pairs
k = 1; %counter for output matrix

for i = 1:pairLength
    if ~strcmp(hamPairs{i,1}, '')  % if not redundant
        outMat{k,1}= hamPairs{i,1};
        outMat{k,2}= hamPairs{i,2};
        outMat{k,3}= hamPairs{i,3};
        k = k+1;
    end
end

%writes to a file to use in cytoscape
XX = fopen([fileName '.txt'],'w');
formatSpec = '%s\t %s\t %d\r\n';
for row = 1:k-1
    if (outMat{row,3}>=hammdist)
    fprintf(XX,formatSpec, outMat{row,:});
    end
end
fclose(XX);
end 

    

  