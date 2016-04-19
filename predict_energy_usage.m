% GP Prediction from PMML file
% Load a GP model from the models folder an then use this model to 
% predict the energy consumption of a certain cut

%prediction part
clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
addpath(genpath('lib'));
load('data/blind_test_accurate.mat')
D = [Data];

% Preload all of the models from PMML
% This is a slow operation so we only want to do it once
% Code block written by Max Ferguson April 14th 2016
for i=1:length(E_blind)
    process_id = ID_blind(i);
    if process_id<0 || process_id >=10
        continue; % Invalid model number
    elseif exist('models','var') && length(models)>=process_id && ~isempty(models{process_id})
        continue % Already loaded this model
    end
    % Load the correct model file from PMML
    filename = sprintf('models/energy-model-%i.pmml',process_id);
    fprintf('Loading model from %s\n',filename);
    models{process_id} = pmml.GaussianProcess(filename);
end


% extract fetures
energy = cell2mat(D(:,9)); %answer
duration = cell2mat(D(:,10)); %duration
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
volume_cut = cell2mat(D(:,29)); %Depth of cut

ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end

%cut_method=zeros(length(y),1);
for i=1:length(feed)

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

%cut_direction=zeros(length(y),1);
for i=1:length(feed)

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

%opeartion 
for i=1:length(feed)
                 
    if strcmp(D{i,33},'Face Milling')
        type_operation(i,1) = 1;      
    elseif strcmp(D{i,33},'Contouring')
        type_operation(i,1) = 2;       
    elseif strcmp(D{i,33},'Slotting')
        type_operation(i,1) = 3;
    elseif strcmp(D{i,33},'Pocketing')
        type_operation(i,1) = 4;
    elseif strcmp(D{i,33},'Spiraling')
        type_operation(i,1) = 5;
    elseif strcmp(D{i,33},'Drilling')
        type_operation(i,1) = 6;             
    else
        type_operation(i,1) = 0;
    end    
end

%label
for i=1:length(feed)

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
    elseif strcmp(D{i,30},'Rapid retract')
        label(i,1) = 8;        
    else
        label(i,1) = 0;
    end   
end


for i=1:length(feed)
    
    % cut with feed
    if (label(i) == 1)
        if (type_operation(i) == 1 ) %face milling
            ID(i,1) = 1;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 2;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 3;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 4;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 5;
        else
            ID(i,1) = 0; %etc
        end  
    elseif (label(i) == 2 && type_operation(i) == 6) %air cut in X-Y 
        ID(i,1) = 6;
    elseif (label(i) == 2 && type_operation(i) ~= 6) %air cut in X-Y 
        ID(i,1) = 7;
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 8;
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 9;        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 9;   
    % aux   
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 10;    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 11;
    else
        ID(i,1) = 12; %no-impact
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
density=energy./length_cut_XYZ;


%total data set
X_blind = input;
ID_blind = input(:,12);
E_blind = energy;
for i=1:length(E_blind)
    if ID(i) == 10 %dwell
        L_blind(i,1)=1;
    else
        L_blind(i,1)=X_blind(i,11);  
    end
end
Y_blind = E_blind./L_blind;
T_blind = input(:,13);


% %clean up data
clean_up_index = find(Y_blind > 0 & Y_blind < inf);
X_blind = X_blind(clean_up_index,:);
ID_blind = ID_blind(clean_up_index);
L_blind = L_blind(clean_up_index);
E_blind = E_blind(clean_up_index);
Y_blind = Y_blind(clean_up_index);
T_blind = T_blind(clean_up_index);





for i=1:length(E_blind)
    % Select the right prediction function
    process_id = ID_blind(i); 

    if (process_id>0 && process_id <10)
        % Select the correct model
        model = models{process_id};
        x_new = X_blind(i,Feature_set{process_id});
        % Generate scores using the PMML model
        [Y_predict(i,1), S_predict(i,1) ] = model.score(x_new);
        E_predict(i,1) = Y_predict(i,1).*L_blind(i);
    else
        Y_predict(i,1)=0;
        S_predict(i,1)=0;
        E_predict(i,1)=0;     
        Y_blind(i,1)=0;
        S_blind(i,1)=0;
        E_blind(i,1)=0;                       
    end
end

number = length(T_blind);
duration = mean(T_blind);
RAE_density = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
predicted_energy = sum(Y_predict.*L_blind);
measured_energy = sum(E_blind);
standard = sqrt(sum(S_predict.*L_blind.^2));
RTE_density = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);

display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])

Y_blind_1 = Y_blind;
E_blind_1 = E_blind;
L_blind_1 = L_blind;

E_predict_1 = E_predict;
S_predict_1 = S_predict;
Y_predict_1 = Y_predict;




























%prediction part
clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
load('data/blind_test_intermediate.mat')
D = [Data];

%extract fetures
energy = cell2mat(D(:,9)); %answer
duration = cell2mat(D(:,10)); %duration
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
volume_cut = cell2mat(D(:,29)); %Depth of cut

ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end

%cut_method=zeros(length(y),1);
for i=1:length(feed)

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

%cut_direction=zeros(length(y),1);
for i=1:length(feed)

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

%opeartion 
for i=1:length(feed)
                 
    if strcmp(D{i,33},'Face Milling')
        type_operation(i,1) = 1;      
    elseif strcmp(D{i,33},'Contouring')
        type_operation(i,1) = 2;       
    elseif strcmp(D{i,33},'Slotting')
        type_operation(i,1) = 3;
    elseif strcmp(D{i,33},'Pocketing')
        type_operation(i,1) = 4;
    elseif strcmp(D{i,33},'Spiraling')
        type_operation(i,1) = 5;
    elseif strcmp(D{i,33},'Drilling')
        type_operation(i,1) = 6;             
    else
        type_operation(i,1) = 0;
    end    
end

%label
for i=1:length(feed)

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
    elseif strcmp(D{i,30},'Rapid retract')
        label(i,1) = 8;        
    else
        label(i,1) = 0;
    end   
end


for i=1:length(feed)
    
    % cut with feed
    if (label(i) == 1)
        if (type_operation(i) == 1 ) %face milling
            ID(i,1) = 1;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 2;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 3;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 4;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 5;
        else
            ID(i,1) = 0; %etc
        end  
    elseif (label(i) == 2 && type_operation(i) == 6) %air cut in X-Y 
        ID(i,1) = 6;
    elseif (label(i) == 2 && type_operation(i) ~= 6) %air cut in X-Y 
        ID(i,1) = 7;
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 8;
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 9;        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 9;   
    % aux   
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 10;    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 11;
    else
        ID(i,1) = 12; %no-impact
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
density=energy./length_cut_XYZ;


%total data set
X_blind = input;
ID_blind = input(:,12);
E_blind = energy;
for i=1:length(E_blind)
    if ID(i) == 10 %dwell
        L_blind(i,1)=1;
    else
        L_blind(i,1)=X_blind(i,11);  
    end
end
Y_blind = E_blind./L_blind;
T_blind = input(:,13);


% %clean up data
clean_up_index = find(Y_blind > 0 & Y_blind < inf);
X_blind = X_blind(clean_up_index,:);
ID_blind = ID_blind(clean_up_index);
L_blind = L_blind(clean_up_index);
E_blind = E_blind(clean_up_index);
Y_blind = Y_blind(clean_up_index);
T_blind = T_blind(clean_up_index);





for i=1:length(E_blind)
    % Select the right prediction function
    process_id = ID_blind(i); 

    if (process_id>0 && process_id <10)
        % Select the correct model file from PMML
        model = models{process_id};
        x_new = X_blind(i,Feature_set{process_id});
        % Generate scores using the PMML model
        [Y_predict(i,1), S_predict(i,1) ] = model.score(x_new);
        E_predict(i,1) = Y_predict(i,1).*L_blind(i);
    else
        Y_predict(i,1)=0;
        S_predict(i,1)=0;
        E_predict(i,1)=0;     
        Y_blind(i,1)=0;
        S_blind(i,1)=0;
        E_blind(i,1)=0;                       
    end
end

number = length(T_blind);
duration = mean(T_blind);
RAE_density = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
predicted_energy = sum(Y_predict.*L_blind);
measured_energy = sum(E_blind);
standard = sqrt(sum(S_predict.*L_blind.^2));
RTE_density = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);

display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])

Y_blind_2 = Y_blind;
E_blind_2 = E_blind;
L_blind_2 = L_blind;

E_predict_2 = E_predict;
S_predict_2 = S_predict;
Y_predict_2 = Y_predict;























%prediction part
clear cut_method cut_direction type_operation label ID input E_predict S_predict Y_predict L_blind Y_blind X_blind
load('data/blind_test_bad.mat')
D = [Data];

%extract fetures
energy = cell2mat(D(:,9)); %answer
duration = cell2mat(D(:,10)); %duration
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
volume_cut = cell2mat(D(:,29)); %Depth of cut

ratio_cut = actual_length_cut./length_cut_XY;
for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end

%cut_method=zeros(length(y),1);
for i=1:length(feed)

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

%cut_direction=zeros(length(y),1);
for i=1:length(feed)

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

%opeartion 
for i=1:length(feed)
                 
    if strcmp(D{i,33},'Face Milling')
        type_operation(i,1) = 1;      
    elseif strcmp(D{i,33},'Contouring')
        type_operation(i,1) = 2;       
    elseif strcmp(D{i,33},'Slotting')
        type_operation(i,1) = 3;
    elseif strcmp(D{i,33},'Pocketing')
        type_operation(i,1) = 4;
    elseif strcmp(D{i,33},'Spiraling')
        type_operation(i,1) = 5;
    elseif strcmp(D{i,33},'Drilling')
        type_operation(i,1) = 6;             
    else
        type_operation(i,1) = 0;
    end    
end

%label
for i=1:length(feed)

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
    elseif strcmp(D{i,30},'Rapid retract')
        label(i,1) = 8;        
    else
        label(i,1) = 0;
    end   
end


for i=1:length(feed)
    
    % cut with feed
    if (label(i) == 1)
        if (type_operation(i) == 1 ) %face milling
            ID(i,1) = 1;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 2;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 3;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 4;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 5;
        else
            ID(i,1) = 0; %etc
        end  
    elseif (label(i) == 2 && type_operation(i) == 6) %air cut in X-Y 
        ID(i,1) = 6;
    elseif (label(i) == 2 && type_operation(i) ~= 6) %air cut in X-Y 
        ID(i,1) = 7;
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 8;
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 9;        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 9;   
    % aux   
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 10;    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 11;
    else
        ID(i,1) = 12; %no-impact
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
density=energy./length_cut_XYZ;


%total data set
X_blind = input;
ID_blind = input(:,12);
E_blind = energy;
for i=1:length(E_blind)
    if ID(i) == 10 %dwell
        L_blind(i,1)=1;
    else
        L_blind(i,1)=X_blind(i,11);  
    end
end
Y_blind = E_blind./L_blind;
T_blind = input(:,13);


% %clean up data
clean_up_index = find(Y_blind > 0 & Y_blind < inf);
X_blind = X_blind(clean_up_index,:);
ID_blind = ID_blind(clean_up_index);
L_blind = L_blind(clean_up_index);
E_blind = E_blind(clean_up_index);
Y_blind = Y_blind(clean_up_index);
T_blind = T_blind(clean_up_index);





for i=1:length(E_blind)
    % Select the right prediction function
    process_id = ID_blind(i); 

    if (process_id>0 && process_id <10)
        % Select the correct model file from PMML
        model = models{process_id};
        x_new = X_blind(i,Feature_set{process_id});
        % Generate scores using the PMML model
        [Y_predict(i,1), S_predict(i,1) ] = model.score(x_new);
        E_predict(i,1) = Y_predict(i,1).*L_blind(i);
    else
        Y_predict(i,1)=0;
        S_predict(i,1)=0;
        E_predict(i,1)=0;     
        Y_blind(i,1)=0;
        S_blind(i,1)=0;
        E_blind(i,1)=0;                       
    end
end

number = length(T_blind);
duration = mean(T_blind);
RAE_density = mean(abs(Y_predict.*L_blind-Y_blind.*L_blind))/mean(Y_blind.*L_blind);
predicted_energy = sum(Y_predict.*L_blind);
measured_energy = sum(E_blind);
standard = sqrt(sum(S_predict.*L_blind.^2));
RTE_density = (sum(Y_predict.*L_blind)-sum(E_blind))/sum(E_blind);

display([number,duration,RAE_density*100,predicted_energy/1000,measured_energy/1000,standard,RTE_density*100])


Y_blind_3 = Y_blind;
E_blind_3 = E_blind;
L_blind_3 = L_blind;

E_predict_3 = E_predict;
S_predict_3 = S_predict;
Y_predict_3 = Y_predict;




figure(1)
hold
stairs(E_predict_1,'b')
stairs(E_blind_1,'r')
xlabel('NC block')
ylabel('Energy $F(x)$','Interpreter','Latex')
legend('Predicted','Measured')
box on
set(gcf,'position',[100,100,500,250])
ylim([0,700])

figure(2)
hold
stairs(E_predict_2,'b')
stairs(E_blind_2,'r')
xlabel('NC block')
ylabel('Energy $F(x)$','Interpreter','Latex')
legend('Predicted','Measured')
box on
set(gcf,'position',[100,100,500,250])
ylim([0,700])

figure(3)
hold
stairs(E_predict_3,'b')
stairs(E_blind_3,'r')
xlabel('NC block')
ylabel('Energy $F(x)$','Interpreter','Latex')
legend('Predicted','Measured')
box on
set(gcf,'position',[100,100,500,250])
ylim([0,700])


Range = [-200:10:200]


figure(4)
hold
hist(E_predict_1-E_blind_1,Range)
xlabel('$|\hat E^i - E^i|$','Interpreter','Latex')
ylabel('Frequency of occurence')
xlim([-200,200])
ylim([0,70])
box on

figure(5)
hold
hist(E_predict_2-E_blind_2,Range)
xlabel('$|\hat E^i - E^i|$','Interpreter','Latex')
ylabel('Frequency of occurence')
xlim([-200,200])
ylim([0,70])
box on

figure(6)
hold
hist(E_predict_3-E_blind_3,Range)
xlabel('$|\hat E^i - E^i|$','Interpreter','Latex')
ylabel('Frequency of occurence')
xlim([-200,200])
ylim([0,70])
box on








%show the total energy prediction
figure(10)
hold
M_1 = sum(E_predict_1);
S_1  = sqrt(sum(S_predict_1.*L_blind_1.^2));
M_2 = sum(E_predict_2);
S_2  = sqrt(sum(S_predict_2.*L_blind_2.^2));
M_3 = sum(E_predict_3);
S_3  = sqrt(sum(S_predict_3.*L_blind_3.^2));


XX_1=linspace(M_1-0.2*M_1,M_1+0.2*M_1,1000);
YY_1 = normpdf(XX_1,M_1,S_1);
plot(XX_1,YY_1,'g','Linewidth',2)
True_1 = normpdf(sum(E_blind_1),M_1,S_1)
plot(sum(E_blind_1), True_1,'o','Markeredgecolor','g','Markerfacecolor','g','Linewidth',2,'Markersize',12')


XX_2=linspace(M_2-0.2*M_2,M_2+0.2*M_2,1000);
YY_2 = normpdf(XX_2,M_2,S_2);
plot(XX_2,YY_2,'b','Linewidth',2)
True_2 = normpdf(sum(E_blind_2),M_2,S_2)
plot(sum(E_blind_2), True_2,'s','Markeredgecolor','b','Markerfacecolor','b','Linewidth',2,'Markersize',12')


XX_3=linspace(M_3-0.2*M_3,M_3+0.2*M_3,1000);
YY_3 = normpdf(XX_3,M_3,S_3);
plot(XX_3,YY_3,'r','Linewidth',2)
True_3 = normpdf(sum(E_blind_3),M_3,S_3)
plot(sum(E_blind_3), True_3,'d','Markeredgecolor','r','Markerfacecolor','r','Linewidth',2,'Markersize',12')


xlabel('Total Energy')
ylabel('PDF')
legend('Prediction (test1)','Measured(test1)','Prediction (test2)','Measured(test2)','Prediction (test3)','Measured(test3)')
box on



line([sum(E_blind_1),sum(E_blind_1)], [0,True_1],'Linestyle','--','color','g','Linewidth',1)
line([sum(E_blind_2),sum(E_blind_2)], [0,True_2],'Linestyle','--','color','b','Linewidth',1)
line([sum(E_blind_3),sum(E_blind_3)], [0,True_3],'Linestyle','--','color','r','Linewidth',1)
xlim([19000,25000])

% 
% 
% RAE_density = mean(abs(Y_predict-Y_test))/mean(Y_test);
% %RAE_density = std(Y_predict-Y_test)/mean(Y_test);
% RAE_energy = mean(abs(E_predict-E_test))/mean(E_test);
% RTE = (sum(E_predict)-sum(E_test))/sum(E_test);
% 
% Estimated_total_energy(1,I) = sum(E_predict);
% True_total_energy(1,I) = sum(E_test);
% Sigma_energy(1,I)  = sqrt((S_predict)'*L_test.^2);
% 







