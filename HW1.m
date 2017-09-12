%Akash Mitra
%am132

% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
x = '3'; y= '5'; %strings
x = 3; y = '5'; %mixed

%Code
if all(ismember(x, '0123456789+-.eEdD'))
    x = str2double(x)
else
    x=x
end
    
if all(ismember(y, '0123456789+-.eEdD'))
    y = str2double(y)
else
    y = y
end

%Answer
x+y


%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length
Nuc = ['A','C','T','G'];
NucPos =randi(4,1,N);
rand_seq = string(Nuc(NucPos)')
rand_seq = strjoin(rand_seq)
rand_seq = strrep(rand_seq,' ','');

% part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

start = strfind(rand_seq, 'ATG');
stop1 = strfind(rand_seq, 'TAA');
stop2 = strfind(rand_seq, 'TGA');
stop3 = strfind(rand_seq, 'TAG');

stop_position = cat(2, stop1, stop2, stop3);
orf2 = []
    for i = 1:length(start)
        orf = min(stop_position(stop_position > start(i)) - start(i))
        orf2 = [orf2; orf]
        longest_length = max(orf2) 
    end
longest_length

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.

N = 500; % define sequence length
Nuc = ['A','C','T','G'];

count = 0;
for i = 1:1000
    NucPos =randi(4,1,N);
    rand_seq = string(Nuc(NucPos)')
    rand_seq = strjoin(rand_seq)
    rand_seq = strrep(rand_seq,' ','');
    start = strfind(rand_seq, 'ATG');
    stop1 = strfind(rand_seq, 'TAA');
    stop2 = strfind(rand_seq, 'TGA');
    stop3 = strfind(rand_seq, 'TAG');

    stop_position = cat(2, stop1, stop2, stop3);
    orf2 = []
    for j = 1:length(start)
        orf = min(stop_position(stop_position > start(j)) - start(j))
        orf2 = [orf2; orf]
    end
    longest_orf = max(orf2);
    if longest_orf >= 50
        count = count+1;
    end
end
prob_greater_50 = count/1000

%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 

mat=[]
N = randi([50,5000],1,100); 
for k = 1:length(N)
    Nuc = ['A','C','T','G'];
    count = 0;
    for i = 1:1000
        NucPos =randi(4,1,N(k));
        rand_seq = string(Nuc(NucPos)');
        rand_seq = strjoin(rand_seq);
        rand_seq = strrep(rand_seq,' ','');
        start = strfind(rand_seq, 'ATG');
        stop1 = strfind(rand_seq, 'TAA');
        stop2 = strfind(rand_seq, 'TGA');
        stop3 = strfind(rand_seq, 'TAG');

        stop_position = cat(2, stop1, stop2, stop3);
        orf2 = []
        for j = 1:length(start)
            orf = min(stop_position(stop_position > start(j)) - start(j));
            orf2 = [orf2; orf];
        end
        longest_orf = max(orf2);
        if longest_orf >= 50 
            count = count+1;
        end
    end
    prob_greater_50 = count/1000
    mat1 = [N(k) prob_greater_50]
    mat=[mat; mat1]
end

scatter(mat(:,1), mat(:,2));
title('Probability of ORF > 50 bp');
xlabel('Length');
ylabel('Probability');

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

% Ans
% When N is large, the ORF would have a higher probability of being greater
% than 50 bp. When N is small, the ORF would have a lower probability of
% being > 50 bp.
% The curve would be sigmoidal. At the lower values of N, the probability
% values would be low. When N gets greater, the probability would be
% tending to 1. At a certain point, the curve would plataeu off at a
% probability of 1, and hence, would appear sigmoidal.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

%{ fid=fopen('qPCRdata.txt', 'r');

fid = fopen('qPCRdata.txt')
head = fgetl(fid);
data = textscan(fid, '%*s %*s %*s %*s %*s %f %*[^\n]','HeaderLines',1);
fid = fclose(fid);
cp_data = cell2mat(data);
l = find(~cp_data)
cp_data = cp_data(1:l-1)

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 

cp_array = reshape(cp_data, 6,12)

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

cp_mean = [] % account for triplicates
for o = 0:3    
   cp_mean(:,o+1) = mean(cp_array(:,o*3+1:(o+1)*3),2); 
end


for m = (1:6) %traverse along conditions
    for n = (1:3) %traverse along genes
        norm_cp(m,n) = 2^(cp_mean(1,n) - cp_mean(m,n) - (cp_mean(1,4) - cp_mean(m,4)));
    end
end

bar(norm_cp);
title("Change in CP between 3 genes for 6 conditions");
xlabel("Conditions");
ylabel("2^ CP");
legend("Gene1","Gene2","Gene3");





%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


