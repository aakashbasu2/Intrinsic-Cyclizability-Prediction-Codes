% This code predicts DNA cyclizability according to the linear physical
% model that consideres contributions from all NN donucleotides.

% Input: Seq (a cell array)
% Output: FF (a matrix)

% Seq is a cell array. Each element is a DNA sequence. All DNA sequences in
% Seq are of the same length. This code predicts the intrinsic
% cyclizability along every sequence and stores it in an array called FF.
% The number of rows in FF is equal to the number of sequences in Seq. For
% each row in FF, i.e. for each sequence in Seq, the columns of FF contain
% the predicted intrinsic cyclizabilities along the sequence. 

%For every sequence, i.e. in every row of FF, the number in the first
%column is the predicted intrinsic cyclizability of the 50 bp DNA fragment
%that spans from position 1 to 50 of the corresponding DNA sequence. The
%second number is the predicted intrinsic cyclizability of the 50 bp DNA
%fragment that spans from position 8 to 57 of the corresponding DNA
%sequence, and so on.

b_linear=[-0.3890
    0.0030
   -0.0173
    0.0414
   -0.0056
    0.0040
    0.0320
    0.0047
    0.0286
   -0.0387
   -0.0185
    0.0305
   -0.0282
    0.0560
   -0.0232
    0.0204];

ctr=[];
LH=cell(1,1);
ctr=1;
flag=1;
counts=floor((numel(Seq{1}) - 50)./7 + 1);
for i=1:numel(Seq)
    t=Seq{i};
    if(numel(t)==923)
        for j=1:7:numel(t)-50
            LH{ctr}=t(j:j+49);
            ctr=ctr+1;
        end
    
    else
       flag=flag+1;
    end
        
end

setsize=numel(LH);
NN=zeros([setsize 16]);
ctr=1;
temp=zeros([1 16]);
for i=1:setsize
    t=LH{i};
    ctr=1;
    for v1=['A','C','G','T']
        for v2=['A','C','G','T']    
            el=[v1,v2];
            temp(ctr)=numel(strfind(t,el));
            ctr=ctr+1;
        end
    end
    NN(i,:)=temp;
end

temp=NN;
temp(:,16)=[];
NNreg=[zeros([setsize 1])+1,temp];
Fp=NNreg*b_linear;
FF=reshape(Fp,[counts,(numel(Seq)-flag+1)]);
    FF=FF';