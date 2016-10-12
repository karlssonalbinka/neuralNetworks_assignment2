% Task 3
clc
clear all

pattern = load('data_task3.txt');
% plot(pattern(:,2), pattern(:,3), 'o')

%%

my = 0.02;
nbrP = length(pattern);
iteration = 10^5;
% iteration = 10^4;
repetitions = 20;
% radialFunc = 5;   %a)
% radialFunc = 20;    %b)
radialFunc = [1:20];  %c)
cError = zeros(1,repetitions);
cErrorFinal = zeros(1,length(radialFunc));

for k_iteration = 1:length(radialFunc)
    k = radialFunc(k_iteration)
    
    g = zeros(k, 1);
    w2_save = zeros(k, repetitions);
    
    %Do everything 20 times & save setup with minimum cError
    for i_time = 1:repetitions
        i_time
        
        w = rand(k,2)*2 - 1;
        t_index = randsample(nbrP, round(nbrP*0.7));
        v_index = (ismember([1:nbrP], t_index) - 1)*-1;
        t_pattern = pattern(t_index,:);
        v_pattern = pattern( logical(v_index),: );
        
        %phase one, competitive learning
        for t = 1:iteration
            index = ceil(rand(1)*nbrP);
            tmp = 0;
            for i = 1:k
                g(i) = exp(-0.5*sum((pattern(index, 2:3) - w(i,:)).^2));
                tmp = tmp + g(i);
            end
            g = g/tmp;
            
            [~, i_minDist] = max(g);
            w(i_minDist,:) = w(i_minDist,:) + my*(pattern(index, 2:3) - w(i_minDist,:));
            
            %Drawing to see that it works
%             if (mod(t, 100)==0)
%                 hold off
%                 figure(1)
%                 plot(pattern(:,2), pattern(:,3), 'o')
%                 hold on
%                 plot(w(:,1), w(:,2), 'r*')
%                 drawnow
%             end
        end
        % ############### Second phase ########################
        % #####################################################
        % Start on the second phase
        %initiate weights and thresholds
        w2 = rand(k,1)*2-1;
        T = rand(1)*2-1;    %threshold
        my2 = 0.1;
        iteration2 = 3000;
        B = 0.5;
        cError_continuous=[];
        g_all = zeros(k, length(t_pattern));
        
        %test
        W = zeros(length(t_pattern), k);
        
        for t = 1:iteration2
            %get random initial pattern
            index = ceil(rand(1)*length(t_pattern));
            
            %calculate g (the input for our system)
            for i = 1:k
                g(i) = exp(-0.5*sum((t_pattern(index, 2:3) - w(i,:)).^2));
            end
            g = g/sum(g);
            
            %get output
            b = w2'*g - T;
            Output = tanh(B*b);
            dT = B*(t_pattern(index,1) - Output)*(1-Output^2);
            dw2 = dT*g;
            T = T - my2*dT;
            w2 = w2 + my2*dw2;
%             W(t,:) = w2';
            
            %calculate classification error continuously (only debugging)
            %         if (mod(t, 50) == 0)
            %
            %             for i = 1:length(t_pattern)
            %                 for j = 1:k
            %                     g_all(j,i) = exp(-0.5*sum((t_pattern(i, 2:3) - w(j,:)).^2));
            %                 end
            %                 g_all(:,i) = g_all(:,i)/sum(g_all(:,i));
            %             end
            %
            %             cError_continuou = [cError_continuou, sum(abs( t_pattern(:,1)' - sign(tanh(B*(w2'*g_all + T))) )) / (2*length(t_pattern))];
            %         end
        end
        % figure(2)   %for seeing cError_continuou and Weights evolution
        % plot(cError_continuou)
        % figure(3)
        % plot(W)
        
        %Calculate cError for each 20 times, save the minimum one
        for i = 1:length(t_pattern)
            for j = 1:k
                g_all(j,i) = exp(-0.5*sum((t_pattern(i, 2:3) - w(j,:)).^2));
            end
            g_all(:,i) = g_all(:,i)/sum(g_all(:,i));
        end
        cError(i_time) = sum(abs( t_pattern(:,1)' - sign(tanh(B*(w2'*g_all + T))) )) / (2*length(t_pattern));
        w2_save(:, i_time) = w2;
        w_save{i_time} = w;
    end
    cErrorFinal(k_iteration) = sum(cError)/length(cError);
end

%% This plots the best sytem
[~, best_i] = min(cError);
best_w2 = w2_save(:,best_i);
best_w = cell2mat(w_save(best_i));

%Draw where the system outputs zero.

nbrPoints = 5*10^3;
x = [25-rand(nbrPoints,1)*40, 15-rand(nbrPoints,1)*25];

disp 1
g_all = zeros(k, length(x));
for i=1:length(x)
    for j = 1:k
        g_all(j,i) = exp(-0.5*sum((x(i, :) - best_w(j,:)).^2));
    end
    g_all(:,i) = g_all(:,i)/sum(g_all(:,i));
end
y = sign(tanh(B*(best_w2'*g_all + T)));

figure(1)
hold on
plot(pattern(:,2), pattern(:,3), 'ob')
plot(x(y==1,1), x(y==1,2), 'g.');
plot(x(y==-1,1), x(y==-1,2), 'y.');
plot(best_w(:,1), best_w(:,2), 'r*')

figure(2)
plot(cErrorFinal);

