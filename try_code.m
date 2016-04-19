clear all
plot_variance = @(x,lower,upper,color) set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor',color);

load('Training1.mat')
load('Training2.mat')
load('Training3.mat')
load('Training4.mat')
load('Training5.mat')
load('Training6.mat')
load('Training7.mat')
load('Training8.mat')
load('Training9.mat')
load('Training10.mat')
load('Training11.mat')
load('Training12.mat')
load('Training13.mat')
load('Training14.mat')
load('Training15.mat')
load('Training16.mat')
load('Training17.mat')
load('Training18.mat')

%D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9];
%D = [Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];
D = [Training1;Training2;Training3;Training4;Training5;Training6;Training7;Training8;Training9;Training10;Training11;Training12;Training13;Training14;Training15;Training16;Training17;Training18];


%D = [Training1;Training2;Training3];
%D = [Training1;Training2;Training3;Training10;Training11;Training12;Training13];

%D=[Training10;Training11;Training12;Training13;Training14;Training15;Training18];



Process_id = 11;
Process_time_cut = 3;



meanfunc = {@meanSum, {@meanLinear, @meanConst}}; %hyp.mean = zeros(size(X,2)+1,1);
likfunc = @likGauss; sn = 1; hyp.lik = log(sn);
covfunc = @covSEard; 
hyp.cov = [log(1000),log(6000),log(1),log(10),log(10),0]';

%1. Cut with feed : Length of cut in X, Length of cut in Y, Volume of material cut, Feed rate, Spindle speed, Depth of Cut, Cutting strategy
%2. Plunge : Length of cut in Z, Volume of material cut, Feed rate, Spindle speed, Depth of Cut
%3. Air-cut : Length of cut in X, Length of cut in Y, Feed rate, Spindle speed
%4. Dwell : Duration, Spindle Speed
%5. Rapid : Length of cut in X, Length of cut in Y, Feed rate, Spindle speed
%6. Other (Auxiliary) : Spindle speed

%% extract the field

% time; %time 
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



ratio_cut = actual_length_cut./length_cut_XY;

for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end



ratio_volum = volume_cut./(length_cut_XY.*depth_cut*9.18);

for i=1:length(ratio_cut)
    if ratio_cut(i)>0;
        ratio_cut(i)=ratio_cut(i);
    else
        ratio_cut(i)=0;
    end
end



%cut_method=zeros(length(y),1);
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


%cut_direction=zeros(length(y),1);
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

% 
%opeartion 
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


%label
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



%prediction function ID

%cut with feed

%face = 11
%contouring = 12
%slotting = 13
%pocketing = 14
%spirating = 15
%plunge = 16




%air cut
%air cut in X Y = 21
%air cut in Z = 22

%aux
%Dwell =31
%Rapid motion = 32
%Others = 33





for i=1:length(energy)
    
    % cut with feed
    if (label(i) == 1) %(cut in x-y direction)

        if (type_operation(i) == 2)     % contouring
            ID(i,1) = 12;
        elseif (type_operation(i) == 3) % splitting
            ID(i,1) = 13;
        elseif (type_operation(i) == 4) % pocketing
            ID(i,1) = 14;
        elseif (type_operation(i) == 5) % spiraling
            ID(i,1) = 15;
        else
            ID(i,1) = 11; % face-milling
        end
        
    elseif (label(i) == 2) %plunge with feed (cut in z direction)
        ID(i,1) = 16;
          
        
        
        
    % air cut      
    elseif (label(i) == 3) % air cut in X-Y
        ID(i,1) = 21;
        
    elseif (label(i) == 4) % air cut in Z
        ID(i,1) = 22;
        
    elseif (label(i) == 5) % air cut in Z
        ID(i,1) = 22;
        
        
    % aux    
    elseif (label(i) == 6) %dwelling
        ID(i,1) = 31;
        
    elseif (label(i) == 7) %rapid motion
        ID(i,1) = 32;
    
    else
        ID(i,1) = 41; %no-labeling
        
    end
    
end
    
  
       

%1: feed
%2: spindle
%3: depth_cut
%4: cut_strategy (cut_method)
%5: volume_cut
%6: direction
%7: length_cut_X (dx)
%8: length_cut_Y (dy)
%9: length_cut_Z (dz)
%10: opearation_type
%11: category

    
          input= [feed,...,
                 spindle,...,
                 depth_cut,...,
                 cut_direction,...,
                 cut_method,...,
                 volume_cut,...,
                 length_cut_X,...,
                 length_cut_Y,...,
                 length_cut_Z,...,
                 length_cut_XY,...,
                 length_cut_XYZ,...,
                 ID,....,
                 duration];
output =energy;
density=energy./length_cut_XY;




%filttering based on job id and duration
index_job = find(ID==Process_id & duration > Process_time_cut  );
[~,random_index] = sort(randn(length(index_job),1));



               %feed spindle depth_cut cut_direction cut_method
index_input = [1 2 3 4 5];
%selection data
X=input(index_job,index_input);
T=input(index_job,13);
L=input(index_job,10);
E=energy(index_job);
Y=E./L;

%shuffling data

X=X(random_index,:);
L=L(random_index,:);
Y=Y(random_index);
E=E(random_index);

%dividing trainig and test data
No_training = round(length(X)*0.9);
X_training = X(1:No_training,:);
Y_training = Y(1:No_training);

X_test = X(No_training+1:end,:);
L_test = L(No_training+1:end,:);
Y_test = Y(No_training+1:end);
E_test = E(No_training+1:end);





%Learning GP
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; %hyp.mean = zeros(size(X,2)+1,1);
likfunc = @likGauss; sn = .1; hyp.lik = log(sn);
covfunc = @covSEard; 
%hyp.cov = zeros(size(X,2)+1,1);

hyp2 = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, X_training, Y_training);

%Predicting GP
[Y_predict,S_predict] = gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, X_test);
E_predict = Y_predict.*L_test;




% %Marginal Uncertainty
% X1 = linspace(mean(X(:,1)),max(X(:,1)),100);
% X2 = linspace(mean(X(:,2)),max(X(:,2)),100);
% X3 = linspace(mean(X(:,3)),max(X(:,3)),100);
% X4 = linspace(mean(X(:,4)),max(X(:,4)),100);
% X5 = linspace(mean(X(:,5)),max(X(:,5)),100);
% 
% YY_predict=zeros(100,100,100,100,100);
% SS_predict=zeros(100,100,100,100,100);

% for i=1:length(X1)
%     for j=1:length(X2)
%         for k=1:length(X3)
%             for p=1:length(X4)
%                 for q=1:length(X5)
%                     [YY_predict(i,j,k,p,q),SS_predict(i,j,k,p,q)] = gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, [X1(i),X2(j),X3(k),X4(p),X5(q)]  );
%                 end
%             end
%         end
%     end
% end

















figure(1)
subplot(2,1,1)
hold
stairs(Y_test,'r')
stairs(Y_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on
subplot(2,1,2)
hold
stairs(E_test,'r')
stairs(E_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on

%show the total energy prediction
figure(2)
hold
M = sum(E_predict);
S  = sqrt(sum(S_predict.*L_test.^2));
XX=linspace(sum(E_predict)-0.1*sum(E_predict),sum(E_predict)+0.1*sum(E_predict),1000)
YY = normpdf(XX,M,S)
plot(XX,YY,'b','Linewidth',2)
XS=linspace(M-1.96*S,M+1.96*S,1000);
YS_upper=normpdf(XS,M,S);
YS_lower=zeros(1,length(XS));
plot_variance(XS,YS_lower,YS_upper,'b')
line([sum(E_test),sum(E_test)],[0,0.001],'color','r','Linewidth',2,'Linestyle','--')
alpha(0.3)
xlabel('Total Energy')
ylabel('PDF')
legend('PDF','95% confidence bound','Measured')

box on





RAE_density = mean(abs(Y_predict-Y_test))/mean(Y_test)
RAE_energy = mean(abs(E_predict-E_test))/mean(E_test)



feed_range = 0:1:1000;
feed = X(:,1);
spindle = X(:,2);
depth_cut = X(:,3);
cut_direction = X(:,4);
cut_method = X(:,5);
density = Y;


%cut == 1          input(:,4)        input(:,5)         input(:,3)    input(:,2)     input(:,2)    input(:,12)
index_11_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle >1400 & spindle <1600 );
index_11_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle >1400 & spindle <1600 );
index_11_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle >1400 & spindle <1600 );
index_11_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle >1400 & spindle <1600 );

index_12_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle >2900 & spindle <3100 );
index_12_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle >2900 & spindle <3100 );
index_12_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle >2900 & spindle <3100 );
index_12_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle >2900 & spindle <3100 );

index_13_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle >4400 & spindle <4600 );
index_13_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle >4400 & spindle <4600 );
index_13_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle >4400 & spindle <4600 );
index_13_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle >4400 & spindle <4600 );

index_14_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle >5900 & spindle <6100 );
index_14_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle >5900 & spindle <6100 );
index_14_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle >5900 & spindle <6100 );
index_14_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle >5900 & spindle <6100 );

         %feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11_1 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_2 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_11_3 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_4 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_12_1 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_2 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_12_3 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_4 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_13_1 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_2 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_13_3 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_4 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_14_1 = [feed_range',6000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_14_2 = [feed_range',6000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_14_3 = [feed_range',6000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_14_4 = [feed_range',6000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];




figure(2)
subplot(3,4,1)
hold
plot(feed(index_11_1),density(index_11_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)


plot(feed(index_11_3),density(index_11_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)

legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])




subplot(3,4,2)
hold
plot(feed(index_12_1),density(index_12_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)


plot(feed(index_12_3),density(index_12_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)

legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(3,4,3)
hold
plot(feed(index_13_1),density(index_13_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)


plot(feed(index_13_3),density(index_13_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)

legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(3,4,4)
hold
plot(feed(index_14_1),density(index_14_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_14_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)


plot(feed(index_14_3),density(index_14_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_14_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)

legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])








%cut == 1.5
index_21_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1.5 & spindle >1400 & spindle <1600 );
index_21_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1.5 & spindle >1400 & spindle <1600 );
index_21_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1.5 & spindle >1400 & spindle <1600 );
index_21_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1.5 & spindle >1400 & spindle <1600 );

index_22_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1.5 & spindle >2900 & spindle <3100 );
index_22_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1.5 & spindle >2900 & spindle <3100 );
index_22_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1.5 & spindle >2900 & spindle <3100 );
index_22_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1.5 & spindle >2900 & spindle <3100 );

index_23_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1.5 & spindle >4400 & spindle <4600 );
index_23_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1.5 & spindle >4400 & spindle <4600 );
index_23_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1.5 & spindle >4400 & spindle <4600 );
index_23_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1.5 & spindle >4400 & spindle <4600 );

index_24_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1.5 & spindle >5900 & spindle <6100 );
index_24_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1.5 & spindle >5900 & spindle <6100 );
index_24_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1.5 & spindle >5900 & spindle <6100 );
index_24_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1.5 & spindle >5900 & spindle <6100 );


         %feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_21_1 = [feed_range',1500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_21_2 = [feed_range',1500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_21_3 = [feed_range',1500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_21_4 = [feed_range',1500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_22_1 = [feed_range',3000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_22_2 = [feed_range',3000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_22_3 = [feed_range',3000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_22_4 = [feed_range',3000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_23_1 = [feed_range',4500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_23_2 = [feed_range',4500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_23_3 = [feed_range',4500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_23_4 = [feed_range',4500*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_24_1 = [feed_range',6000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_24_2 = [feed_range',6000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_24_3 = [feed_range',6000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_24_4 = [feed_range',6000*ones(length(feed_range),1), 1.5*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];




subplot(3,4,5)
hold
plot(feed(index_21_1),density(index_21_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_21_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_21_3),density(index_21_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_21_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(3,4,6)
hold
plot(feed(index_22_1),density(index_22_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_22_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_22_3),density(index_22_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_22_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(3,4,7)
hold
plot(feed(index_23_1),density(index_23_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_23_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_23_3),density(index_23_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_23_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(3,4,8)
hold
plot(feed(index_24_1),density(index_24_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_24_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_24_3),density(index_24_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_24_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])








%cut == 3
index_31_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 3 & spindle >1400 & spindle <1600 );
index_31_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 3 & spindle >1400 & spindle <1600 );
index_31_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 3 & spindle >1400 & spindle <1600 );
index_31_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 3 & spindle >1400 & spindle <1600 );

index_32_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 3 & spindle >2900 & spindle <3100 );
index_32_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 3 & spindle >2900 & spindle <3100 );
index_32_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 3 & spindle >2900 & spindle <3100 );
index_32_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 3 & spindle >2900 & spindle <3100 );

index_33_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 3 & spindle >4400 & spindle <4600 );
index_33_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 3 & spindle >4400 & spindle <4600 );
index_33_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 3 & spindle >4400 & spindle <4600 );
index_33_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 3 & spindle >4400 & spindle <4600 );

index_34_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 3 & spindle >5900 & spindle <6100 );
index_34_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 3 & spindle >5900 & spindle <6100 );
index_34_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 3 & spindle >5900 & spindle <6100 );
index_34_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 3 & spindle >5900 & spindle <6100 );


XX_31_1 = [feed_range',1500*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_31_2 = [feed_range',1500*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_31_3 = [feed_range',1500*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_31_4 = [feed_range',1500*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_32_1 = [feed_range',3000*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_32_2 = [feed_range',3000*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_32_3 = [feed_range',3000*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_32_4 = [feed_range',3000*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_33_1 = [feed_range',4500*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_33_2 = [feed_range',4500*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_33_3 = [feed_range',4500*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_33_4 = [feed_range',4500*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_34_1 = [feed_range',6000*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_34_2 = [feed_range',6000*ones(length(feed_range),1), 3*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_34_3 = [feed_range',6000*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_34_4 = [feed_range',6000*ones(length(feed_range),1), 3*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];




subplot(3,4,9)
hold
plot(feed(index_31_1),density(index_31_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_31_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_31_3),density(index_31_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_31_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])

subplot(3,4,10)
hold
plot(feed(index_32_1),density(index_32_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_32_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_32_3),density(index_32_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_32_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])


subplot(3,4,11)
hold
plot(feed(index_33_1),density(index_33_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_33_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_33_3),density(index_33_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_33_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])


subplot(3,4,12)
hold
plot(feed(index_34_1),density(index_34_1),'go')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_34_1)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','g')
alpha(0.3)
plot(feed(index_34_3),density(index_34_3),'rd')
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_34_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('X-cut(data)','X-cut(predict)','Y-cut(data)','Y-cut(predict)')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])













































%Test for blind test
clear cut_method cut_direction type_operation label ID input
load('blind_test_accurate.mat')
%load('blind_test_intermediate.mat')
%load('blind_test_bad.mat')
D = [Data];


%% extract the field

% time; %time 
energy = cell2mat(D(:,9)); %answer


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

% 
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
        if (type_operation(i) == 1) %face milling
            ID(i,1) = 11;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 12;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 13;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 14;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 15;
        else
            ID(i,1) = 17; %etc
        end
        
    elseif (label(i) == 2) %plunge with feed (cut in z direction)
        ID(i,1) = 16;
          
        
        
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 21;
        
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 22;
        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 22;
        
        
    % aux    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 31;
        
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 32;
    
    else
        ID(i,1) = 41; %no-impact
        
    end
    
end
    
  
input= [feed,...,
     spindle,...,
     depth_cut,...,
     cut_direction,...,
     cut_method,...,
     volume_cut,...,
     length_cut_X,...,
     length_cut_Y,...,
     length_cut_Z,...,
     length_cut_XY,...,
     length_cut_XYZ,...,
     ID];



output =energy;
density=energy./length_cut_XYZ;




%face milling
index_job = find(input(:,12)==Process_id);


               %feed spindle depth_cut cut_direction cut_method
index_input = [1 2 3 4 5];
%selection data
X_blind=input(index_job,index_input);
L_blind=input(index_job,11);
E_blind=energy(index_job);
Y_blind=E_blind./L_blind;

[Y_blind_predict,S_blind_predict] = gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind);
E_blind_predict = Y_blind_predict.*L_blind;



% figure(3)
% subplot(2,2,1)
% hold
% stairs(Y_blind,'r')
% stairs(Y_blind_predict,'b')
% legend('True','Predicted')
% subplot(2,2,3)
% stairs(Y_blind_predict-Y_blind,'k')
% legend('Error')
% subplot(2,2,2)
% hold
% stairs(E_blind,'r')
% stairs(E_blind_predict,'b')
% legend('True','Predicted')
% subplot(2,2,4)
% stairs(E_blind_predict-E_blind,'k')
% legend('Error')
% 
% RAE_blind_density = mean(abs(Y_blind_predict-Y_blind))/mean(Y_blind)
% RAE_blind_energy = mean(abs(E_blind_predict-E_blind))/mean(E_blind)





figure(4)
subplot(2,1,1)
hold
stairs(Y_blind,'r')
stairs(Y_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on
subplot(2,1,2)
hold
stairs(E_blind,'r')
stairs(E_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on

%show the total energy prediction
figure(5)
hold
M = sum(E_blind_predict);
S  = sqrt(sum(S_blind_predict.*L_blind.^2));
XX=linspace(sum(E_blind_predict)-0.2*sum(E_blind_predict),sum(E_blind_predict)+0.2*sum(E_blind_predict),1000)
YY = normpdf(XX,M,S)
plot(XX,YY,'b','Linewidth',2)
XS=linspace(M-1.96*S,M+1.96*S,1000);
YS_upper=normpdf(XS,M,S);
YS_lower=zeros(1,length(XS));
plot_variance(XS,YS_lower,YS_upper,'b')
line([sum(E_blind),sum(E_blind)],[0,0.001],'color','r','Linewidth',2,'Linestyle','--')
alpha(0.3)
xlabel('Total Energy')
ylabel('PDF')
legend('PDF','95% confidence bound','Measured')

box on



figure(6)
feed_range = 0:1:1000;

%cut == 1
index_11_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 1450 & spindle < 1550 & ID == Process_id);
index_11_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 1450 & spindle < 1550 & ID == Process_id);
index_11_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 1450 & spindle < 1550 & ID == Process_id);
index_11_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 1450 & spindle < 1550 & ID == Process_id);

index_12_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 2950 & spindle < 3050 & ID == Process_id);
index_12_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 2950 & spindle < 3050 & ID == Process_id);
index_12_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 2950 & spindle < 3050 & ID == Process_id);
index_12_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 2950 & spindle < 3050 & ID == Process_id);

index_13_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 4450 & spindle < 4550 & ID == Process_id);
index_13_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 4450 & spindle < 4550 & ID == Process_id);
index_13_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 4450 & spindle < 4550 & ID == Process_id);
index_13_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 4450 & spindle < 4550 & ID == Process_id);



         %feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11_1 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_2 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_11_3 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_4 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_12_1 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_2 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_12_3 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_4 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_13_1 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_2 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_13_3 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_4 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];





subplot(2,3,1)
hold
plot(feed(index_11_3),density(index_11_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,2)
hold
plot(feed(index_12_3),density(index_12_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,3)
hold
plot(feed(index_13_3),density(index_13_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,4)
hold
plot(feed(index_11_4),density(index_11_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,5)
hold
plot(feed(index_12_4),density(index_12_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,6)
hold
plot(feed(index_13_4),density(index_13_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])

























%Test for blind test
clear cut_method cut_direction type_operation label ID input
load('blind_test_intermediate.mat')

D = [Data];


%% extract the field

% time; %time 
energy = cell2mat(D(:,9)); %answer


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

% 
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
        if (type_operation(i) == 1) %face milling
            ID(i,1) = 11;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 12;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 13;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 14;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 15;
        else
            ID(i,1) = 17; %etc
        end
        
    elseif (label(i) == 2) %plunge with feed (cut in z direction)
        ID(i,1) = 16;
          
        
        
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 21;
        
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 22;
        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 22;
        
        
    % aux    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 31;
        
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 32;
    
    else
        ID(i,1) = 41; %no-impact
        
    end
    
end
    
  
input= [feed,...,
     spindle,...,
     depth_cut,...,
     cut_direction,...,
     cut_method,...,
     volume_cut,...,
     length_cut_X,...,
     length_cut_Y,...,
     length_cut_Z,...,
     length_cut_XY,...,
     length_cut_XYZ,...,
     ID];



output =energy;
density=energy./length_cut_XYZ;




%face milling
index_job = find(input(:,12)==Process_id);


               %feed spindle depth_cut cut_direction cut_method
index_input = [1 2 3 4 5];
%selection data
X_blind=input(index_job,index_input);
L_blind=input(index_job,11);
E_blind=energy(index_job);
Y_blind=E_blind./L_blind;

[Y_blind_predict,S_blind_predict] = gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind);
E_blind_predict = Y_blind_predict.*L_blind;



figure(7)
subplot(2,1,1)
hold
stairs(Y_blind,'r')
stairs(Y_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on
subplot(2,1,2)
hold
stairs(E_blind,'r')
stairs(E_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on

%show the total energy prediction
figure(8)
hold
M = sum(E_blind_predict);
S  = sqrt(sum(S_blind_predict.*L_blind.^2));
XX=linspace(sum(E_blind_predict)-0.2*sum(E_blind_predict),sum(E_blind_predict)+0.2*sum(E_blind_predict),1000)
YY = normpdf(XX,M,S)
plot(XX,YY,'b','Linewidth',2)
XS=linspace(M-1.96*S,M+1.96*S,1000);
YS_upper=normpdf(XS,M,S);
YS_lower=zeros(1,length(XS));
plot_variance(XS,YS_lower,YS_upper,'b')
line([sum(E_blind),sum(E_blind)],[0,0.001],'color','r','Linewidth',2,'Linestyle','--')
alpha(0.3)
xlabel('Total Energy')
ylabel('PDF')
legend('PDF','95% confidence bound','Measured')

box on




%Blind test intermediat
%spindle = 1700, 2800, 4300

figure(9)
feed_range = 0:1:1000;

%cut == 1
index_11_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 1650 & spindle < 1750 & ID == Process_id);
index_11_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 1650 & spindle < 1750 & ID == Process_id);
index_11_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 1650 & spindle < 1750 & ID == Process_id);
index_11_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 1650 & spindle < 1750 & ID == Process_id);

index_12_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 2750 & spindle < 3850 & ID == Process_id);
index_12_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 2750 & spindle < 3850 & ID == Process_id);
index_12_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 2750 & spindle < 3850 & ID == Process_id);
index_12_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 2750 & spindle < 3850 & ID == Process_id);

index_13_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 4250 & spindle < 4350 & ID == Process_id);
index_13_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 4250 & spindle < 4350 & ID == Process_id);
index_13_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 4250 & spindle < 4350 & ID == Process_id);
index_13_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 4250 & spindle < 4350 & ID == Process_id);



         %feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11_1 = [feed_range',1700*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_2 = [feed_range',1700*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_11_3 = [feed_range',1700*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_4 = [feed_range',1700*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_12_1 = [feed_range',2800*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_2 = [feed_range',2800*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_12_3 = [feed_range',2800*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_4 = [feed_range',2800*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_13_1 = [feed_range',4300*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_2 = [feed_range',4300*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_13_3 = [feed_range',4300*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_4 = [feed_range',4300*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];





subplot(2,3,1)
hold
plot(feed(index_11_3),density(index_11_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,2)
hold
plot(feed(index_12_3),density(index_12_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,3)
hold
plot(feed(index_13_3),density(index_13_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,4)
hold
plot(feed(index_11_4),density(index_11_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,5)
hold
plot(feed(index_12_4),density(index_12_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,6)
hold
plot(feed(index_13_4),density(index_13_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])






























%Blind test bad
%spindle = 2130, 2400, 3750
%Test for blind test
clear cut_method cut_direction type_operation label ID input
load('blind_test_bad.mat')
D = [Data];


%% extract the field

% time; %time 
energy = cell2mat(D(:,9)); %answer


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

% 
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
        if (type_operation(i) == 1) %face milling
            ID(i,1) = 11;
        elseif (type_operation(i) == 2) %contouring
            ID(i,1) = 12;
        elseif (type_operation(i) == 3) %splitting
            ID(i,1) = 13;
        elseif (type_operation(i) == 4) %pocketing
            ID(i,1) = 14;
        elseif (type_operation(i) == 5) %spiralling
            ID(i,1) = 15;
        else
            ID(i,1) = 17; %etc
        end
        
    elseif (label(i) == 2) %plunge with feed (cut in z direction)
        ID(i,1) = 16;
          
        
        
    % air cut      
    elseif (label(i) == 3) %air cut in X-Y 
        ID(i,1) = 21;
        
    elseif (label(i) == 4) %air cut in Z
        ID(i,1) = 22;
        
    elseif (label(i) == 5) %air cut in Z
        ID(i,1) = 22;
        
        
    % aux    
    elseif (label(i) == 6) %dwell
        ID(i,1) = 31;
        
    elseif (label(i) == 7 ||label(i) == 8) %rapid motion
        ID(i,1) = 32;
    
    else
        ID(i,1) = 41; %no-impact
        
    end
    
end
    
  
input= [feed,...,
     spindle,...,
     depth_cut,...,
     cut_direction,...,
     cut_method,...,
     volume_cut,...,
     length_cut_X,...,
     length_cut_Y,...,
     length_cut_Z,...,
     length_cut_XY,...,
     length_cut_XYZ,...,
     ID];



output =energy;
density=energy./length_cut_XYZ;




%face milling
index_job = find(input(:,12)==Process_id);


               %feed spindle depth_cut cut_direction cut_method
index_input = [1 2 3 4 5];
%selection data
X_blind=input(index_job,index_input);
L_blind=input(index_job,11);
E_blind=energy(index_job);
Y_blind=E_blind./L_blind;

[Y_blind_predict,S_blind_predict] = gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, X_blind);
E_blind_predict = Y_blind_predict.*L_blind;



% figure(3)
% subplot(2,2,1)
% hold
% stairs(Y_blind,'r')
% stairs(Y_blind_predict,'b')
% legend('True','Predicted')
% subplot(2,2,3)
% stairs(Y_blind_predict-Y_blind,'k')
% legend('Error')
% subplot(2,2,2)
% hold
% stairs(E_blind,'r')
% stairs(E_blind_predict,'b')
% legend('True','Predicted')
% subplot(2,2,4)
% stairs(E_blind_predict-E_blind,'k')
% legend('Error')
% 
% RAE_blind_density = mean(abs(Y_blind_predict-Y_blind))/mean(Y_blind)
% RAE_blind_energy = mean(abs(E_blind_predict-E_blind))/mean(E_blind)





figure(10)
subplot(2,1,1)
hold
stairs(Y_blind,'r')
stairs(Y_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on
subplot(2,1,2)
hold
stairs(E_blind,'r')
stairs(E_blind_predict,'b')
xlabel('NC block')
ylabel('Energy per length')
legend('Measured','Predicted')
box on

%show the total energy prediction
figure(11)
hold
M = sum(E_blind_predict);
S  = sqrt(sum(S_blind_predict.*L_blind.^2));
XX=linspace(sum(E_blind_predict)-0.2*sum(E_blind_predict),sum(E_blind_predict)+0.2*sum(E_blind_predict),1000)
YY = normpdf(XX,M,S)
plot(XX,YY,'b','Linewidth',2)
XS=linspace(M-1.96*S,M+1.96*S,1000);
YS_upper=normpdf(XS,M,S);
YS_lower=zeros(1,length(XS));
plot_variance(XS,YS_lower,YS_upper,'b')
line([sum(E_blind),sum(E_blind)],[0,0.001],'color','r','Linewidth',2,'Linestyle','--')
alpha(0.3)
xlabel('Total Energy')
ylabel('PDF')
legend('PDF','95% confidence bound','Measured')

box on




%spindle = 2130, 2400, 3750
figure(12)
feed_range = 0:1:1000;

%cut == 1
index_11_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 2080 & spindle < 2180 & ID == Process_id);
index_11_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 2080 & spindle < 2180 & ID == Process_id);
index_11_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 2080 & spindle < 2180 & ID == Process_id);
index_11_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 2080 & spindle < 2180 & ID == Process_id);

index_12_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 2350 & spindle < 2450 & ID == Process_id);
index_12_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 2350 & spindle < 2450 & ID == Process_id);
index_12_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 2350 & spindle < 2450 & ID == Process_id);
index_12_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 2350 & spindle < 2450 & ID == Process_id);

index_13_1 = find(cut_direction==1 & cut_method==1 & depth_cut == 1 & spindle > 3700 & spindle < 3800 & ID == Process_id);
index_13_2 = find(cut_direction==1 & cut_method==2 & depth_cut == 1 & spindle > 3700 & spindle < 3800 & ID == Process_id);
index_13_3 = find(cut_direction==2 & cut_method==1 & depth_cut == 1 & spindle > 3700 & spindle < 3800 & ID == Process_id);
index_13_4 = find(cut_direction==2 & cut_method==2 & depth_cut == 1 & spindle > 3700 & spindle < 3800 & ID == Process_id);



         %feed            spindle                          depth_cut                     cut_direction                   cut_method
XX_11_1 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_2 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_11_3 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_11_4 = [feed_range',1500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_12_1 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_2 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_12_3 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_12_4 = [feed_range',3000*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];

XX_13_1 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_2 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),1*ones(length(feed_range),1), 2*ones(length(feed_range),1)];
XX_13_3 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 1*ones(length(feed_range),1)];
XX_13_4 = [feed_range',4500*ones(length(feed_range),1), 1*ones(length(feed_range),1),2*ones(length(feed_range),1), 2*ones(length(feed_range),1)];





subplot(2,3,1)
hold
plot(feed(index_11_3),density(index_11_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,2)
hold
plot(feed(index_12_3),density(index_12_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,3)
hold
plot(feed(index_13_3),density(index_13_3),'bs', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_3)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','b')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,4)
hold
plot(feed(index_11_4),density(index_11_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_11_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,5)
hold
plot(feed(index_12_4),density(index_12_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_12_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])



subplot(2,3,6)
hold
plot(feed(index_13_4),density(index_13_4),'rd', 'Markersize',10, 'Linewidth',2)
[m,s]=gp(hyp2, @infExact, [], covfunc, likfunc, X_training, Y_training, XX_13_4)
plot_variance(feed_range,(m-sqrt(s))',(m+sqrt(s))','r')
alpha(0.3)
legend('Measurement','Prediction')
xlabel('Feed rate')
ylabel('Energy/length cut')
xlim([0,700])
ylim([0,15])


