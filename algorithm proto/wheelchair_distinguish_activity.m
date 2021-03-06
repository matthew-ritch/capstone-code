clc
clear
close all

datafolder="sensor_tests_activities";
%%categories={"inactive", "general", "lift", "good", "loop", "choppy"};
categories={'inactive', "general", 'lift', 'good', 'choppy'};
%categories={"inactive", "general", "lift", "choppy", "good"};


%% read file names
data_files = dir(strcat(datafolder, "/sensors_recording*"));
cal_file = dir(strcat(datafolder, "/*calibration*"));

%% get calibration file
[~,~,calib]=xlsread(strcat(cal_file.folder, "\", cal_file.name));
%parse calibration file
[cal_acc,cal_mag,cal_fs]=parseData(calib);
cal_acc=mean(cal_acc);
%set proper magnitude but use direction
%cal_acc=9.81*cal_acc/norm(cal_acc);
cal_mag=mean(cal_mag);
cal_mag=cal_mag/norm(cal_mag);

samples ={};
labels=[];

if ~exist('data.mat')|| false
    for i =1:length(data_files)
        for j=1:length(categories)
            type(j)=~isempty(strfind(data_files(i).name,categories{j}));
        end
        type=find(type);
        %read in data file
        [~,~,data]=xlsread(strcat(data_files(i).folder, "\", data_files(i).name));
        %drop first and last few seconds
        data=data(2000:end-2000,:);
        %parse
        [acc,mag,fs]=parseData(data);
        len=256/fs; %seconds for each sample
        
        length_t=(data{end,1}-data{1,1})/(10^9);
        t=0:1/fs:(1/fs)*(size(acc,1)-1);
        %tell user file details
        fprintf("file: %s, fs: %d, length: %d\n\n", data_files(i).name, fs, length_t)    
        %split data into samples of time-length len seconds
        dataArray = splitData (acc,mag,fs,len,t, cal_acc, cal_mag);

        samples=[samples,dataArray];
        labels=[labels, type*ones(1,length(dataArray))];
    end
    %save file so don't have to read every time
    save('data.mat','samples','labels')
else
    %load file
    a=load('data.mat');
    samples=a.samples;
    labels=a.labels;
end


features=[];
for i=1:length(samples)
    %isolate one sample of length len seconds
    trial = samples(:,i);
    %make starting matrix of accelerometer and magnetometer readings by axis
    matrix = [trial{1}(:,1), trial{1}(:,2), trial{1}(:,3),trial{2}(:,1), trial{2}(:,2), trial{2}(:,3)];
    %make 0-centered matrix without drift for integral and derivative calculations
    matrix_centered = matrix - mean(matrix);
    %take derivative and integral
    derivs=100*diff(matrix_centered);
    derivs=[0,0,0,0,0,0;derivs];
    ints = (1/100)*cumtrapz(matrix_centered);
    %concatenate these readings
    trip=[matrix, derivs, ints];
    %get correlation coefficients of these columns
    R = corrcoef (trip);
    R = R(triu(ones(size(R,1),size(R,2)),1)==1)';
    %only use acceleration (x,y,z) correlations
    R = R(1:(48));
    %add to full feature vector
    features = [features;mean(trip), rms(trip), std(trip), R];
    %features = [features;mean(trip),  std(trip)];
end
labels=labels';
lab_back=labels;
feat_back=features;

%% do simple active vs inactive
cats=categories(1:2);
cats{1}='Not Actively Pushing';
cats{2}='Actively Pushing';
features=feat_back;
labels=lab_back;
labels(labels~=1)=2;
rng(111)
%randomize order
rv=randperm(length(labels));
features=features(rv,:);
labels=labels(rv);
%split dataset into test and train 60%-40%
rVec=rand(1,size(features,1));
training=rVec<.6;
testing=rVec>=.6;
X=features(training,:);
X_labs=labels(training);
Y=features(testing,:);
Y_labs=labels(testing);
X(:,end+1)=1;
Y(:,end+1)=1;
%do linear regression
B = X\X_labs;
X_pred=X*B;
correct=round(X_pred)==X_labs;
trainaccuracyReg_simp=sum(correct)/length(correct);
%apply params to check accuracy
Y_pred=Y*B;
correct=round(Y_pred)==Y_labs;
testaccuracyReg_simp=sum(correct)/length(correct);
%plot results
figure()
scatter(Y_labs(correct), Y_pred(correct),'k');
hold on
scatter(Y_labs(~correct), Y_pred(~correct),'r');
xticks([1,2])
xlim([0,3])
xticklabels({"Not Actively Pushing","Actively Pushing"})
yticks([1,2])
ylim([0,3])
yticklabels({"Not Actively Pushing","Actively Pushing"})
title("Detecting User Activity")
n=length(Y_labs);
ps=[round(100*sum(Y_labs==1)/n),round(100*sum(Y_labs~=1)/n)];
accuracy = round(100*sum(correct)/n);
n
fprintf("Detecting User Activity Results:\nData breakdown: %i%% inactive and %i%% active\n%i%% Accuracy.\n\n\n",ps(1),ps(2),accuracy)
B=B';
for i = 1:length(B)
    fprintf(upper(sprintf("%d, ",B(i))))
end
fprintf('\n')

figure()
Y_labs=cats(Y_labs);
Y_pred=(round(Y_pred));
Y_pred(Y_pred==0)=1;
Y_pred=cats(Y_pred);

cm=confusionchart(Y_labs, Y_pred);
title("Detecting User Activity")



%% discriminate types of activity
cats=categories(3:end);
features=feat_back;
labels=lab_back;
%take only activity types
features=features((labels~=1)&(labels~=2),:);
labels=labels((labels~=1)&(labels~=2),:);
labels=labels-2; %make labels start at 1 (we removed first two categories)
%randomize order
rng(111)
rv=randperm(length(labels));
features=features(rv,:);
labels=labels(rv);
%split dataset into testing/training
rVec=rand(1,size(features,1));
training=rVec<.6;
testing=rVec>=.6;
X=features(training,:);
X_labs=labels(training);
Y=features(testing,:);
Y_labs=labels(testing);
X(:,end+1)=1;
Y(:,end+1)=1;
%do linear regression
B = X\X_labs;
X_pred=X*B;
correct=round(X_pred)==X_labs;
trainaccuracyReg_active=sum(correct)/length(correct);
%use regression to predict on testing set
Y_pred=Y*B;
correct=round(Y_pred)==Y_labs;
testaccuracyReg_active=sum(correct)/length(correct);
%plot results
figure()
scatter(Y_labs(correct), Y_pred(correct),'k');
hold on
scatter(Y_labs(~correct), Y_pred(~correct),'r');
hold off
xticks(1:max(Y_labs))
xlim([0,max(Y_labs)+1])
xticklabels({"Lifting", "Good Form",  "Choppy"}) %"Good Form (loop)",
yticks(1:max(Y_labs))
ylim([0,max(Y_labs)+1])
yticklabels({"Lifting", "Good Form",  "Choppy"}) %"Good Form (loop)",
title("Detecting User Activity Type")
n=length(Y_labs);
ns=[sum(Y_labs==1),sum(Y_labs==2),sum(Y_labs==3)];
ps=round(100*ns./n);
n
fprintf("Detecting User Activity Results:\nData breakdown: %i%% Lifting, %i%% Good, and %i%% choppy.\n%i%% Accuracy.\n\n\n",ps(1),ps(2),ps(3),round(100*sum(correct)/n))
B=B';
for i = 1:length(B)
    fprintf(upper(sprintf("%d, ",B(i))))
end
fprintf('\n')
Y_labs=cats(Y_labs);
Y_pred=(round(Y_pred));
Y_pred(Y_pred<1)=1;
Y_pred(Y_pred>3)=3;
Y_pred=cats(Y_pred);

figure()
cm=confusionchart(Y_labs, (Y_pred));
title("Classifying User Activity")



%   
% %% random forest
% t = templateTree('NumVariablesToSample','all',...
%     'PredictorSelection','interaction-curvature','Surrogate','on');
% Mdl = fitcensemble(X,X_labs,'Method','Bag','NumLearningCycles',200, ...
%     'Learners',t);
% 
% yHat = oobPredict(Mdl);
% R2 = corr(Mdl.Y,yHat)^2;
% 
% 
% 
% preds= predict(Mdl, Y);
% 
% correct = preds == Y_labs;
% accuracy = sum(correct)/length(correct);




%%calculate features: https://pdfs.semanticscholar.org/77aa/e609ab8a80f2056cf66588ad5fce31a1410a.pdf
    %try: average mag of acceleration, clockwise vs counterclockwise, stdev, skewness, kurtosis, pairwise correlation of acc, rms values, sum of values 
%features=make_features(samples);

function features=make_features(samples)

    for i=1:length(samples)
        trial = samples{i};
        matrix = [trial{1}(:,1), trial{1}(:,2), trial{1}(:,3),trial{2}(:,1), trial{2}(:,2), trial{2}(:,3)];
    end


end




function dataArray = splitData (acc,mag,fs,len,t, cal_acc, cal_mag)
    dataArray={};
    number = floor((size(acc,1)/fs)/len); %how many samples of length len exist
    [acc, mag] = remove_gravity(acc, mag, cal_acc, cal_mag);
    
    for i=1:number
       start = (i-1)*len;
       fin = (i)*len;
       Index= (t>=start) & (t<fin);
       accThis=acc(Index,:);
       magThis=mag(Index,:);
       dataArray=[dataArray,[{accThis};{magThis}]];
    end
    
end

function [acc, mag] = remove_gravity(acc, mag, cal_acc, cal_mag)
    
    magNorms=sqrt(mag(:,1).^2+ mag(:,2).^2+ mag(:,3).^2);
    mag=mag./magNorms;
    
    %get translation to rotate gravity into new reference orientation for
    %each time point
    %https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    u = mag(:,:);                      % input
    base = cal_mag;                      % original orientation
    R=[];
    I=[1,0,0;0,1,0;0,0,1];
    ssc=@(x)[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ]; %skewsymmetric
    RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);  %rotation calc function
    for j=1:length(u)
       a=u(j,:); %get mag for this time point
       R{j}=RU(a,base); %transposing takes inverse rotation       
    end
    %subtract rotated gravity
    for j=1:length(acc)
        %acc(j,:)=acc(j,:)-(R{j}*cal_acc')';
    end  
end




function [acc,mag,fs] = parseData (data)
    %dump the first 100 (~.2 seconds) and last second
    data=data(100:end,:);
    mags=data(strcmp(data(:,2),'MAG'),:);
    accs=data(strcmp(data(:,2),'ACC'),:);
    mags=cell2mat(mags(:,[1,3,4,5]));
    accs=cell2mat(accs(:,[1,3,4,5]));
    nMags=size(mags,1);
    nAccs=size(accs,1);
    modparam = round(nAccs/nMags); 
    modparam = 5;
    accs=accs(1:modparam:end,:);
    %correct sizes
    while (size(accs,1)>size(mags,1))
        accs=accs(1:end-1,:);
    end
    while (size(mags,1)>size(accs,1))
        mags=mags(1:end-1,:);
    end
    fs=size(accs,1)/((accs(end,1)-accs(1,1))/(10^9));
    acc=accs(:,[2,3,4]);
    mag=mags(:,[2,3,4]);
    %plot(diff(mags(:,1)))

end
