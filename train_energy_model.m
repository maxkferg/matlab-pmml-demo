% Create GP models for energy prediction and save them to PMML
% Models are used to predict the energy consumption of a certain cut
%
% Written by Jinkyoo Park, 2015
% Modified by Max Ferguson, 2016

% Add helper functions to the path
addpath(genpath('lib'));

plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color);
process_ID=[1 2 3 4 5 6 7 8 9 10];

% Base covariance function scale parameters
Hyp = [log(1000),log(6000),log(6),log(3),log(3),log(1)];

Feature_set{1} = [1 2 3 4 5 ];
Feature_set{2} = [1 2 3 4 ];
Feature_set{3} = [1 2 3 4 ];
Feature_set{4} = [1 2 3 4 ];
Feature_set{5} = [1 3 4];
Feature_set{6} = [1 2 3 4];
Feature_set{7} = [1 2 3 4];
Feature_set{8} = [1 2 4];
Feature_set{9} = [1];
Feature_set{10} = [2 4];


load('data/Training1.mat')
load('data/Training2.mat')
load('data/Training3.mat')
load('data/Training4.mat')
load('data/Training5.mat')
load('data/Training6.mat')
load('data/Training7.mat')
load('data/Training8.mat')
load('data/Training9.mat')
load('data/Training10.mat')
load('data/Training11.mat')
load('data/Training12.mat')
load('data/Training13.mat')
load('data/Training14.mat')
load('data/Training15.mat')
load('data/Training16.mat')
load('data/Training17.mat')
load('data/Training18.mat')

D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9;Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];


%% Extract the fields
energy = cell2mat(D(:,9)); %energy consumption
duration = cell2mat(D(:,10)); %duration of operation
feed = cell2mat(D(:,11)); %duration of operation
spindle = cell2mat(D(:,12)); %spindle speed
length_cut_X = abs(cell2mat(D(:,19))); %code dx
length_cut_Y = abs(cell2mat(D(:,20))); %code dy
length_cut_Z = abs(cell2mat(D(:,23))); %code dy
length_cut_XY = cell2mat(D(:,21)); %code length_cut
length_cut_XYZ = cell2mat(D(:,31)); %code length_cut
actual_dx = cell2mat(D(:,25)); %actual dx 
actual_dy = cell2mat(D(:,26)); %actual dy
actual_length_cut = sqrt(actual_dx.^2+actual_dy.^2);
depth_cut = cell2mat(D(:,27)); %Depth of cut
area_cut = cell2mat(D(:,28)); %Depth of cut
volume_cut = cell2mat(D(:,29)); %Depth of cut

% Ratio cut
ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end


for i=1:length(energy)

    if strcmp(D{i,24},'Conventional')
        cut_method(i,1) = 1;   
    elseif strcmp(D{i,24},'Climb')
        cut_method(i,1) = 2;
    elseif strcmp(D{i,24},'Both')
        cut_method(i,1)=3;
    else
        cut_method(i,1)=0;
    end 
end


for i=1:length(energy)
    if abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))==0
        cut_direction(i,1) = 1;   
    elseif abs(cell2mat(D(i,19)))==0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,1) = 2;
    elseif abs(cell2mat(D(i,19)))>0 & abs(cell2mat(D(i,20)))>0
        cut_direction(i,1)=3;
    else
        cut_direction(i,1)=0;
    end
        
end

% Operation 
for i=1:length(energy)
         
    if strcmp(D{i,32},'Face Milling')
        type_operation(i,1) = 1;      
    elseif strcmp(D{i,32},'Contouring')
        type_operation(i,1) = 2;       
    elseif strcmp(D{i,32},'Slotting')
        type_operation(i,1) = 3;
    elseif strcmp(D{i,32},'Pocketing')
        type_operation(i,1) = 4;
    elseif strcmp(D{i,32},'Spiraling')
        type_operation(i,1) = 5;
    elseif strcmp(D{i,32},'Drilling')
        type_operation(i,1) = 6;             
    else
        type_operation(i,1) = 0;
    end
        
end

% Label
for i=1:length(energy)
    if strcmp(D{i,30},'Cut with Feed')
        label(i,1) = 1;
    elseif strcmp(D{i,30},'Plunge with feed')
        label(i,1) = 2;
    elseif strcmp(D{i,30},'Air-Cut')
        label(i,1) = 3;
    elseif strcmp(D{i,30},'Air-Cut in Z while plunging')
        label(i,1) = 4;
    elseif strcmp(D{i,30},'Air-cut in Z while retracting')
        label(i,1) = 5;
    elseif strcmp(D{i,30},'Dwell')
        label(i,1) = 6;  
    elseif strcmp(D{i,30},'No Cut - Rapid motion')
        label(i,1) = 7;
    else
        label(i,1) = 0;
    end   
end

for i=1:length(energy)
    % Cut with feed
    if (label(i) == 1) %(cut in x-y direction)
     
        if (type_operation(i) == 1)     % face milling
            ID(i,1) = 1;
        elseif (type_operation(i) == 2) % contouring
            ID(i,1) = 2;
        elseif (type_operation(i) == 3) % splitting
            ID(i,1) = 3;
        elseif (type_operation(i) == 4) % pocketing
            ID(i,1) = 4;
        elseif (type_operation(i) == 5) % spiraling
            ID(i,1) = 5;
        else
            ID(i,1) = 0; % 
        end  
    elseif (label(i)==2 && type_operation(i)==6) %other plunge with feee disregard
        ID(i,1) = 6; %   Drilling
    elseif (label(i)==2 && type_operation(i)~=6) %other plunge with feee disregard
        ID(i,1) = 7; %   Plunge
    % air cut      
    elseif (label(i) == 3) % air cut in X-Y
        ID(i,1) = 8;       
    elseif (label(i) == 4) % air cut in Z
        ID(i,1) = 9;     
    % aux    
    elseif (label(i) == 7) %rapid motion
        ID(i,1) = 10;
    elseif (label(i) == 6) %dwell
        ID(i,1) = 11;
    else
        ID(i,1) = 11; %no-labeling
    end
end
    
input= [feed,...,
        spindle,...,
        depth_cut,...,
        cut_direction,...,
        cut_method,...,
        ratio_cut,...,
        length_cut_X,...,
        length_cut_Y,...,
        length_cut_Z,...,
        length_cut_XY,...,
        length_cut_XYZ,...,
        ID,...,
        duration];
output =energy;

% Total data set
X = input;
E = energy;
for i=1:length(E)
    if ID(i) == 10 %dwell
        L(i,1)=1;
    else
        L(i,1)=X(i,11);  
    end
end
Y = E./L;


% Clean up data
clean_up_index = find(Y > 0 & Y < inf);
X = X(clean_up_index,:);
ID = ID(clean_up_index);
L = L(clean_up_index);
E = E(clean_up_index);
Y = Y(clean_up_index);


% Training time filtering
time_filter_index = find(X(:,13) > 3);
X = X(time_filter_index,:);
ID = ID(time_filter_index);
L = L(time_filter_index);
E = E(time_filter_index);
Y = Y(time_filter_index);


% Training individul prediction function
for I=1:length(process_ID)
    
    feature_index = Feature_set{I};
    process_id = process_ID(I);
    
    %Select data depending on the process ID
    if process_id == 10 %dwell
        index_job{I} = find(X(:,12)==process_id );
    else                                             %time           %non-zero cut
        index_job{I} = find(X(:,12)==process_id & X(:,13)>3 & X(:,11)>0); %cutting related
    end

    % Selection data
    X_training=X(index_job{I},feature_index);
    E_training=E(index_job{I});
    L_training=L(index_job{I});
    T_training=X(index_job{I},13);
    Y_training=E_training./L_training;

    % Setting GP initial parameters 
    hyp.mean = [];
    hyp.cov = [Hyp(feature_index),0];
    meanfunc = 'MeanZero';
    covfunc = 'ARDSquaredExponentialKernel';
    likfunc = 'Gaussian'; sn = 0.1; hyp.lik = log(sn);
    inffunc = 'Exact';
    filename = sprintf('models/energy-model-%i.pmml',I);
    
    % Learning GP model and saving to PMML
    fprintf('Optimizing Gaussian Process model %i\n',I);
    fprintf('Using features %s in model %i\n\n',strjoin(feature_index,','),I);
    p = pmml.GaussianProcess(hyp, inffunc, meanfunc, covfunc, likfunc, X_training, Y_training);
    p.optimize(-100);
    p.toPMML(filename);
    fprintf('Saved Gaussian Process model to %s\n',filename);
end