% This code predicts DNA cyclizability according to the linear physical
% model that consideres contributions from all 16 NN dinucleotide pairs and
% all 136 NN-NN Helica Repetition Indices.

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


b_combined=[ -0.4181
    0.0022
   -0.0160
    0.0439
   -0.0059
    0.0045
    0.0321
    0.0056
    0.0266
   -0.0396
   -0.0190
    0.0310
   -0.0310
    0.0593
   -0.0209
    0.0237
   -0.4086
    0.1999
    0.5685
    1.2294
   -0.2147
    0.2200
    0.2630
    0.1358
    0.1704
    0.3163
   -0.3237
    0.2169
   -0.2394
   -0.6057
   -0.0061
   -0.2298
    0.1055
   -0.0847
    0.4119
   -0.0738
   -0.0435
    0.0762
    0.0524
   -0.0437
   -0.0305
   -0.1000
   -0.0356
   -0.2013
   -0.1948
   -0.0784
    0.0038
    0.3099
    0.6038
   -0.0833
    0.0000
   -0.0756
   -0.1526
   -0.0493
   -0.0967
   -0.2424
   -0.0541
   -0.1634
   -0.3184
   -0.0665
   -0.2026
   -0.7013
   -0.1657
    0.1585
    0.2825
    0.3795
    0.2375
    0.1576
   -0.2146
    0.2689
   -0.1992
   -0.5864
    0.1635
   -0.1493
    0.0655
    0.0819
    0.0245
    0.0116
   -0.1140
    0.1849
    0.0749
   -0.1431
    0.3371
    0.3834
    0.1539
   -0.0851
    0.1194
    0.1648
   -0.0247
   -0.0593
   -0.0811
   -0.2342
   -0.2200
    0.2497
    0.0088
   -0.0302
   -0.2680
   -0.7649
   -0.1661
    0.1816
    0.2788
   -0.0382
   -0.0525
   -0.0503
    0.0196
    0.0017
    0.0273
   -0.6239
    0.0788
    0.0271
    0.0722
   -0.0206
   -0.1957
    0.0965
    0.0135
    0.1795
   -0.4235
   -0.0814
    0.1892
    0.1330
   -0.0063
    0.0375
    0.0772
   -0.1331
   -0.5011
   -0.0624
   -0.1488
    0.0668
    0.1956
    0.0209
   -0.1500
    0.1695
    0.1844
   -0.0942
    0.5073
    0.1604
    0.3890
    0.3239
   -0.3301
    0.0134
   -0.1251
    0.3461
    0.6580
    0.3546
    0.0211
   -0.1454
    0.4560
    0.3443
    0.4055
   -0.5746
    0.2363
    0.1506];


ctr=[];
LH=cell(1,1);
ctr=1;
flag=1;
counts=floor((numel(Seq{1}) - 50)./7 + 1);
for i=1:numel(Seq)
    t=Seq{i};
    if(numel(t)>0)%==2001)
        for j=1:7:numel(t)-49 %%%%%%%%%%
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
    i
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

temp=NN;
temp(:,16)=[];
NNreg=[zeros([setsize 1])+1,temp,metric];
Fp=NNreg*b_combined;
FF=reshape(Fp,[counts,(numel(Seq)-flag+1)]);
FF=FF';