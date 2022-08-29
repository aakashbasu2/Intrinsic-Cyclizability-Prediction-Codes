% This code predicts DNA cyclizability according to the linear physical
% model that consideres contributions from all NN-NN Helical Repetition
% Indices.

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

b_periodic = [-0.1203
   -0.3364
    0.1967
    0.6972
    1.2084
   -0.2556
    0.2450
    0.2020
    0.1286
    0.0986
    0.2706
   -0.3064
    0.2813
   -0.2167
   -0.5983
   -0.0914
   -0.1927
    0.1478
   -0.0244
    0.4155
   -0.1165
   -0.0245
    0.0542
    0.0074
   -0.0896
   -0.0216
   -0.1522
   -0.0489
   -0.1719
   -0.1830
   -0.1581
    0.0077
   -0.3743
    0.6981
    0.0065
    0.1067
   -0.0353
   -0.0344
    0.0897
    0.0207
   -0.2253
    0.0329
   -0.0075
   -0.2033
   -0.0436
   -0.0811
   -0.4802
   -0.2165
    0.2030
    0.2107
    0.3040
    0.1836
    0.1384
   -0.2990
    0.2598
   -0.1644
   -0.5777
    0.0846
   -0.1444
    0.2448
    0.1116
   -0.0431
   -0.0012
   -0.1782
    0.1211
    0.0492
   -0.1121
    0.3591
    0.3485
    0.1149
   -0.1454
   -0.2803
    0.1687
   -0.0228
   -0.0340
   -0.0821
   -0.2028
   -0.1104
    0.4110
    0.0710
   -0.0728
   -0.1405
   -0.3710
   -0.2559
    0.1206
    0.2543
   -0.1298
   -0.0649
   -0.1055
   -0.0267
   -0.1437
    0.0124
   -0.3947
    0.0379
   -0.0266
    0.0202
   -0.0376
   -0.2142
    0.1034
   -0.1117
    0.2118
   -0.2127
   -0.1869
    0.1278
    0.1127
   -0.0228
   -0.0468
    0.0486
   -0.1581
   -0.4093
   -0.0718
   -0.0609
    0.0604
    0.1844
   -0.0253
   -0.1518
    0.4161
    0.2032
   -0.0824
    0.5001
    0.0866
    0.3625
    0.0035
   -0.2204
    0.0270
   -0.1925
    0.5027
    0.2971
    0.3647
   -0.0284
   -0.0904
    0.5350
    0.2979
    0.3731
   -0.1031
    0.1980
    0.0215];
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


k=1;
S=cell(1,1);
for i=['A','T','G','C']
    for j=['A','T','G','C']
        S{k}=[i,j];
        k=k+1;
    end
end


S={'AA','AT','TA','TT','AC','AG','GA','CA','TC','TG','GT','CT','GG','GC','CG','CC'};
mm=zeros([1 18]);
setsize=numel(LH);
mmm=zeros([setsize 6]);
metric=zeros([setsize 136]);
goo=1;
arr=zeros([1 49]);
for n1=1:16
    for n2=n1:16
        pair{goo}=[S{n1},S{n2}];
        AA=[];
        for i=1:setsize
            t=LH{i};
            s1=strfind(t,S{n1});
            s2=strfind(t,S{n2});
            arr(:)=0;
            arr(s1)=1;
            arr(s2)=1;
            ctr=1;
            for r=[9,10,11,19,20,21,29,30,31,4,5,6,14,15,16,24,25,26]
                s=0;
                for j=1:50-r-1
                    s=s+arr(j).*arr(j+r);
                end
                s=s./(50-r-1);
                mm(ctr)=s;
                ctr=ctr+1;
            end
            mmm(i,:)=[max(mm(1:3)),max(mm(4:6)),max(mm(7:9)),max(mm(10:12)),max(mm(13:15)),max(mm(16:18))];
        end
        nummmm=mmm(:,1)+mmm(:,2)+mmm(:,3)-mmm(:,4)-mmm(:,5)-mmm(:,6);
        metric(:,goo)=nummmm;
        goo=goo+1;
    end
end

NNreg=[zeros(setsize,1)+1,metric];
Fp=NNreg*b_periodic;
FF=reshape(Fp,[counts,(numel(Seq)-flag+1)]);
FF=FF';