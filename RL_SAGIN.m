clear all; close all; clc;
tic %啟動計時器

coverage_length=100;   %km
coverage_width=100;    %km
N_Leos=5;   %The Number of Low Erath Orbit
N_Haps=10;  %The Number of Haps
N_GroundUsers=15; %The Number of GroundUsers
Height_Leo=300;  % 300km
Height_Haps=20;  % 20km
Height_Ground=0; % 0 km

Leo_coordinates = [rand(1, N_Leos) * coverage_length; rand(1, N_Leos) * coverage_width; ones(1, N_Leos) * Height_Leo];
Haps_coordinates = [rand(1, N_Haps) * coverage_length; rand(1, N_Haps) * coverage_width; ones(1, N_Haps) * Height_Haps];
GroundUser_coordinates = [rand(1, N_GroundUsers) * coverage_length; rand(1, N_GroundUsers) * coverage_width; ones(1, N_GroundUsers) * Height_Ground];


i = 3; 
x_Leo_i = Leo_coordinates(1, i);
y_Leo_i = Leo_coordinates(2, i);
z_Leo_i = Leo_coordinates(3, i);


frequency = 2.4e9; % 2.4 GHz
c = 3e8; % Speed of light
FSPL_Leo_Haps = zeros(N_Leos, N_Haps);  
FSPL_Haps_Ground = zeros(N_Haps, N_GroundUsers);  
for i = 1:N_Leos
    for j = 1:N_Haps
        distance_Leo_Haps = norm(Leo_coordinates(:, i) - Haps_coordinates(:, j));
        FSPL_Leo_Haps(i, j) = 20*log10(distance_Leo_Haps) + 20*log10(frequency) + 20*log10(4*pi/c) + 147.55;
    end
end
for i = 1:N_Haps
    for j = 1:N_GroundUsers
        distance_Haps_Ground = norm(Haps_coordinates(:, i) - GroundUser_coordinates(:, j));
        FSPL_Haps_Ground(i, j) = 20*log10(distance_Haps_Ground) + 20*log10(frequency) + 20*log10(4*pi/c)+147.55;
    end
end


% Calculate final channel gains
mean_exp = 1; % Mean of the exponential distribution
channel_Leo_Haps = FSPL_Leo_Haps .* exprnd(mean_exp, size(FSPL_Leo_Haps));
channel_Haps_Ground = FSPL_Haps_Ground .* exprnd(mean_exp, size(FSPL_Haps_Ground));


% % LEO to HAPS connection
% LEO_Haps_connection = zeros(N_Leos, N_Haps);   % LEO to HAPS connection matrix
% for i = 1:N_Haps
%     connected_leo = randperm(N_Leos, 1);
%     LEO_Haps_connection(connected_leo, i) = 1;  % Connect the selected LEO to this HAPS
% end
% % HAPS to GroundUser connection
% Haps_GroundUser_connection = zeros(N_Haps, N_GroundUsers);   % HAPS to GroundUser connection matrix
% for i = 1:N_GroundUsers
%     connected_haps = randperm(N_Haps, 1);
%     Haps_GroundUser_connection(connected_haps, i) = 1;  % Connect the selected HAPS to this GroundUser
% end



pts = 10;      %pts可以理解為不同傳輸功率的數量
P_max = 10;     %最大功率選擇點的大小
P = linspace(1,P_max,pts);  %產生等距的數列 (1~P_max) (pts個)
state_level = 10;

R_Leo = 0;
iter = 1e6;
LEO_Haps_connection = zeros(N_Leos, N_Haps);
for i=1:iter
    atmp_leos = ceil(rand(1,N_Leos)*pts);  %產生 1*N_Haps 個 1~pts的整數
    Power_Leos = zeros(1,N_Leos);   %產生 1*N_Leo 個0， 每個設備的發射功率
    for k=1:N_Leos
         Power_Leos(k) = P(atmp_leos(k));
    end
    % LEO to HAPS connection
    LEO_Haps_connection(:)=0;
    for j = 1:N_Haps
        connected_leo = randperm(N_Leos, 1);
        LEO_Haps_connection(connected_leo,j) = 1;  % Connect the selected LEO to this HAPS
    end
    
    Ttmp_leo = reward_v2( Power_Leos, N_Leos,N_Haps,channel_Leo_Haps,LEO_Haps_connection);
    if Ttmp_leo>R_Leo
        R_Leo = Ttmp_leo;
    end
    progress = i / iter;  
    if mod(i, 4000) == 0 || i == iter
    fprintf('Progress (Find R_max_leo): %.2f%%\n', progress * 100);
    end
end
R_Leo_max = R_Leo;

R_Haps = 0;
Haps_GroundUser_connection = zeros(N_Haps, N_GroundUsers);   % HAPS to GroundUser connection matrix
for i=1:iter
    atmp = ceil(rand(1,N_Haps)*pts);  
    Power = zeros(1,N_Haps);  
    for k=1:N_Haps
        Power(k) = P(atmp(k));
    end

    Haps_GroundUser_connection(:)=0;
    for j = 1:N_GroundUsers
        connected_haps = randperm(N_Haps, 1);
        Haps_GroundUser_connection(connected_haps, j) = 1;  % Connect the selected HAPS to this GroundUser
    end
    Ttmp = reward_v2(Power,N_Haps,N_GroundUsers,channel_Haps_Ground,Haps_GroundUser_connection);
    if Ttmp>R_Haps
        R_Haps = Ttmp;
    end
    progress = i / iter;
    if mod(i, 4000) == 0 || i == iter
    fprintf('Progress(Find R_max_Haps): %.2f%%\n', progress * 100);
    end
end
R_Haps_max = R_Haps;


disp("R_LEO_MAX:")
disp(R_Leo_max)
disp("R_HAPS_MAX:")
disp(R_Haps_max)

% 
% R_Leo_max=11;
% R_Haps_max=7;



% Define Q-tables for LEOs and HAPS
Q1_link = cell(1, N_Haps);
for i = 1:N_Haps
    Q1_link{i} = zeros(state_level,N_Leos, N_Leos);
end

Q2_link = cell(1, N_GroundUsers);
for i = 1:N_GroundUsers
    Q2_link{i} = zeros(state_level,N_Haps, N_Haps);
end
% Define epsilon-greedy exploration parameters
Q1_power = cell(1,N_Leos);
for i=1:N_Leos
    Q1_power{i} = zeros(pts,state_level);
end
Q2_power = cell(1,N_Haps);
for i=1:N_Haps
    Q2_power{i} = zeros(pts,state_level);
end

iter=1e6;
eta = 0.05; % learning rate
gamma = 0.3; % next important

s_Leo_power = ones(1,iter);
s_Haps_power= ones(1,iter);
act_Leo_power = ones(N_Leos,iter);
act_Haps_power= ones(N_Haps,iter);
%init
Power_Leos = ones(1, N_Leos) * 10;
Power_Haps = ones(1, N_Haps) * 10;
% Define initial states
a_Haps_link = ones(1, N_Haps);
s_Haps_link = ones(1, N_Haps);
s_Haps_next_link=ones(1, N_Haps);
a_GroundUsers_link = ones(1, N_GroundUsers);
s_GroundUsers_link = ones(1, N_GroundUsers);
s_GroundUsers_next_link = ones(1, N_GroundUsers);

state_level=10;
state = @(x, Rmax) ceil(x*state_level/Rmax);


throughput_Leo= ones(1,iter);
throughput_Haps= ones(1,iter);

for i = 1:iter
    % Select actions
    if rand() < 0.5 || i <=iter/4
        % Random action
        if mod(i,2)==0   %Power
            a_leo = ceil(rand(1,N_Leos)*pts);   
            a_haps= ceil(rand(1,N_Haps)*pts);
            Power_Leo = zeros(1,N_Leos);
            Power_Haps= zeros(1,N_Haps);
            for k=1:N_Leos
                Power_Leo(k) = P(a_leo(k));
            end
            for k=1:N_Haps
                Power_Haps(k)=P(a_haps(k));
            end
        else  %Link
            for j = 1:N_Haps
                a_Haps_link(j) = randi(N_Leos);  % Randomly select a Leo for each Haps
                s_Haps_next_link(j)=a_Haps_link(j);
            end
            for j = 1:N_GroundUsers
                a_GroundUsers_link(j) = randi(N_Haps);  % Randomly select a Haps for each Ground User
                s_GroundUsers_next_link(j)=a_GroundUsers_link(j);
            end
        end
    else % Optimal action based on Q-values
        if mod(i,2)==0   %power
            a_leo = zeros(1,N_Leos);
            a_haps = zeros(1,N_Haps);
            for k=1:N_Leos
                [val ,a_leo(k)] = max(Q1_power{k}(:,s_Leo_power(i)));
            end
            for k=1:N_Haps
                [val ,a_haps(k)] = max(Q2_power{k}(:,s_Haps_power(i)));
            end
            Power_Leo = zeros(1,N_Leos);
            Power_Haps= zeros(1,N_Haps);
            for k=1:N_Leos
                Power_Leo(k) = P(a_leo(k));
            end
            for k=1:N_Haps
                Power_Haps(k)=P(a_haps(k));
            end
        else             %Link
            LEO_Haps_connection = zeros(N_Leos, N_Haps);
            Haps_GroundUser_connection = zeros(N_Haps, N_GroundUsers); 
            for j = 1:N_Haps
                [~, a_Haps_link(j)] = max(Q1_link{j}(s_Leo_power(i),s_Haps_link(j), :));
                s_Haps_next_link(j)=a_Haps_link(j);  % Connect the Haps to the selected Leo
                LEO_Haps_connection(a_Haps_link(j), j) = 1;
            end
            for j = 1:N_GroundUsers
                [~,a_GroundUsers_link(j)] = max(Q2_link{j}(s_Haps_power(i),s_GroundUsers_link(j), :));
                s_GroundUsers_next_link(j)=a_GroundUsers_link(j);  % Connect the Ground User to the selected Haps
                Haps_GroundUser_connection(a_GroundUsers_link(j), j) = 1;
            end
        end
    end

    throughput_Leo(i) = reward_v2(Power_Leos, N_Leos, N_Haps, channel_Leo_Haps, LEO_Haps_connection);
    throughput_Haps(i) = reward_v2(Power_Haps, N_Haps, N_GroundUsers, channel_Haps_Ground, Haps_GroundUser_connection);
    
   	s_Leo_power(i+1) = state(throughput_Leo(i),R_Leo_max);
    s_Haps_power(i+1)= state(throughput_Haps(i),R_Haps_max);
    if s_Leo_power(i+1)<=0
    s_Leo_power(i+1) = 1;
    elseif s_Leo_power(i+1)>state_level
        s_Leo_power(i+1) = state_level;
    end
    if s_Haps_power(i+1)<=0
        s_Haps_power(i+1) = 1;
    elseif s_Haps_power(i+1)>state_level
        s_Haps_power(i+1) = state_level;
    end

    % Update Q-table
    if mod(i,2)==0
        for k=1:N_Leos  
            Q1_power{k}(a_leo(k), s_Leo_power(i)) = Q1_power{k}(a_leo(k), s_Leo_power(i)) + eta*(throughput_Leo(i) + gamma*max(Q1_power{k}(:, s_Leo_power(i+1))) - Q1_power{k}(a_leo(k), s_Leo_power(i)));
        end
        for k=1:N_Haps
            Q2_power{k}(a_haps(k), s_Haps_power(i)) = Q2_power{k}(a_haps(k), s_Haps_power(i)) + eta*(throughput_Haps(i) + gamma*max(Q2_power{k}(:, s_Haps_power(i+1))) - Q2_power{k}(a_haps(k), s_Haps_power(i)));
        end
    else
        for j = 1:N_Haps
            Q1_link{j}(s_Leo_power(i),s_Haps_link(j), a_Haps_link(j)) = Q1_link{j}(s_Leo_power(i),s_Haps_link(j), a_Haps_link(j)) + eta * (throughput_Leo(i) + gamma * max(Q1_link{j}(s_Leo_power(i+1),s_Haps_next_link(j),:)) - Q1_link{j}(s_Leo_power(i),s_Haps_link(j), a_Haps_link(j)));
        end
        for j = 1:N_GroundUsers
            Q2_link{j}(s_Haps_power(i),s_GroundUsers_link(j), a_GroundUsers_link(j)) = Q2_link{j}(s_Haps_power(i),s_GroundUsers_link(j), a_GroundUsers_link(j)) + eta * (throughput_Haps(i) + gamma * max(Q2_link{j}(s_Haps_power(i+1),s_GroundUsers_next_link(j), :)) - Q2_link{j}(s_Haps_power(i),s_GroundUsers_link(j), a_GroundUsers_link(j)));
        end
        s_Haps_link=s_Haps_next_link;
        s_GroundUsers_link=s_GroundUsers_next_link;
    end
    
   
    progress = i / iter  ;
    if mod(i, 4000) == 0 || i == iter
        fprintf('Progress: %.2f%%\n', progress * 100);
    end
end

% Plot Leos information
figure; 
batch = (iter)/1000;
TT = ones(1,(iter)/batch);
for i=1:(iter)/batch
    TT(i) = mean(throughput_Leo((i-1)*batch+1:i*batch));
end
subplot(1,1,1); plot(TT); xlabel('Epoch'); ylabel('Averaged throughput of each 1e3 Epochs')
set(gca,'FontSize',14)
set(gcf,'color','w')

% Plot HAPS information
figure; 
batch = (iter)/1000;
TT_Haps = ones(1,(iter)/batch);
for i=1:(iter)/batch
    TT_Haps(i) = mean(throughput_Haps((i-1)*batch+1:i*batch));
end
subplot(1,1,1); plot(TT_Haps); xlabel('Epoch'); ylabel('Averaged throughput of each 1e3 Epochs')
set(gca,'FontSize',14)
set(gcf,'color','w')

toc %停止計時器
