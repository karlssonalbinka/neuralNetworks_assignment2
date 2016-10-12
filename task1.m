clc
clear all
nbrP = 1000;
nbrW = 100;
dist = zeros(1, nbrW);
dist_w = zeros(1, nbrW);
% rho_0 = 100;
rho_0 = 5;
rho_c=0.9;
my_0 = 0.1;
my_c = 0.01;
tau = 200;
t_o = 10^3;
t_c = 5*10^4;
n_func = @(bestIndex, index, rho) exp(-sum( (bestIndex-index)^2 )/(2*rho^2)); %neighbourhood function


%produce 1000 points in triangle
y0 = sqrt(3);
pattern = [];
while (length(pattern) < nbrP)
   point = rand(1,2);
   if (point(2) < point(1)*y0 && point(2)< y0*(1-point(1)))
       pattern = [pattern; point];
   end
end


% hold on
% plot([0; 0.5], [0; y0*0.5], 'r', [0.5; 1], [y0*0.5; 0], 'r', [0; 1], [0; 0], 'r');
% plot(pattern(:,1), pattern(:,2), '.b');
% axis equal

%%
close all

%initialize weights
w = [];
while (length(w) < nbrW) %initialize weights in same triangle
   point = rand(1,2);
   if (point(2) < point(1)*y0 && point(2)< y0*(1-point(1)))
       w = [w; point];
   end
end

mult = 1; %test scale
pattern=pattern*mult;
w = w*mult;

for n_it = 1:t_o
% for n_it = 1:5
    my = my_0*exp(-n_it/tau);
    rho = rho_0*exp(-n_it/tau);
    
    %draw random pattern and find closest w
    index = ceil(rand(1)*nbrP);
    for i = 1:nbrW
        dist(i) = sqrt( sum((pattern(index,:)-w(i,:)).^2) );
    end
    [~, min_i] = min(dist);
    
    %update weights
    for i = 1:nbrW
%         tmp = n_func(w(i,:), w(min_i,:), rho);
        tmp = n_func(min_i, i, rho);
        w(i,:) = w(i,:) + my*tmp*(pattern(index,:) - w(i,:));
    end
    
    if (mod(n_it, 1) == 0)
        hold off
        plot([0; 0.5]*mult, [0; y0*0.5]*mult, 'r', [0.5; 1]*mult, [y0*0.5; 0]*mult, 'r', [0; 1]*mult, [0; 0], 'r');
        hold on
        plot(pattern(:,1), pattern(:,2), '.g');
        plot(w(:,1), w(:,2), 'b.');
        axis equal
        title('order');
        drawnow
    end
end
%%
% #################################3 CONVERGANCE #########################
my = my_c;
rho = rho_c;
for n_it = 1:t_c

    
    %draw random pattern and find closest w
    index = ceil(rand(1,1)*nbrP);
    for i = 1:nbrW
        dist(i) = sqrt((pattern(index, 1)-w(i,1))^2 + (pattern(index, 2)-w(i,2))^2 );
    end
    [~, min_i] = min(dist);
%     plot(pattern(index,1), pattern(index,2), 'k*');
%     plot(w(min_i,1), w(min_i,2), 'y*');
    
    %update weights
    for i = 1:nbrW
%         w(min_i,:)
%         w(i,:)
        tmp = n_func(min_i, i, rho);
        w(i,:) = w(i,:) + my*tmp*(pattern(index,:) - w(i,:) );
    end
    
    if (mod(n_it, 10) == 0)
        hold off
        plot([0; 0.5]*mult, [0; y0*0.5]*mult, 'r', [0.5; 1]*mult, [y0*0.5; 0]*mult, 'r', [0; 1]*mult, [0; 0], 'r');
        hold on
        plot(pattern(:,1), pattern(:,2), '.g');
        plot(w(:,1), w(:,2), 'b');
        axis equal
        title('convergance');
        drawnow
    end
end



