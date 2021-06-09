%% load data
P = table2array(readtable('Personal_attribute.csv'));
load Trajectory_location.mat
load Trajectory_building_id.mat
buildings = table2array(readtable('Building_attribute.csv'));

%% parameters initialization
simu_weeks = 25;          % simulation weeks
init_infe_rate = 0.001;   % initial infetion rate
latent_period = 2;        % latent period in day
alpha = 0.8;              % vaccine efficacy
v = 0.25;                 % vaccination rate 
s_c = 0.25;               % student composition (% of distance-learning)

% simulation period everday
start_t = 6*30+30; % 6:30
end_t = 23*30+30; % 23:30

% average number of close contacts in each activity environment
near_num_dining_hall = 10;
near_num_lecture_hall = 10;
near_num_library = 5;

% parameters of libaray and dining_hall
libaray_floors = 5; % number of floors in the library
libaray_capa_each = 300; % capacity of each floor in the library
capacity_dining_hall = [1800,600,300,1800,100,200,400]; % capacity of each dining_hall

% infection rate in each activity environment
p_dining_hall = 3.0316e-04;
p_lecture = 3.3037e-05;
p_residential = 1.7375e-04;
p_road = 2.3200e-05;
p_com = 9.5000e-08;
K = 1/3; %Decay coefficient for the infection rate of an Ip-status agent

%% model initialization
% number of students on campus
num = size(P,1)*(1-s_c); 

% extract id of each activity environment
dining_hall_id = buildings(buildings(:,2)==1,1);
lecture_hall_id = buildings(buildings(:,2)==2,:);
residential_id = buildings(buildings(:,2)==3,:);
library_id = buildings(buildings(:,2)==4,:);

% initialize capacity and adjacency matrix in lecture halls
tab = tabulate(P(:,4));
count = tab(:,2);
capacity_lecture = ceil(count/50)*50;
class_p = cell(1,length(capacity_lecture));
class_near = cell(1,length(capacity_lecture));

% initialize adjacency matrix in dining halls
dining_hall_p = cell(1,length(capacity_dining_hall));
dining_hall_near = cell(1,length(capacity_dining_hall));

% initialize adjacency matrix in library
libaray_p = cell(1,libaray_floors);
libaray_near = cell(1,libaray_floors);
    
% remove agents of distance-learning
list = round(1:1/s_c:size(P,1));
P(list,:) = [];
Timetable(list,:,:) = [];
timetable(list,:) = [];
    
% vacination record
temp = randperm(num);
P(temp(1:round(alpha*v*num)),5) = 1;
    
% introduce initial infected agents randomly 
num_list = 1:num;
num_list(P(:,5) == 1)=[]; % exclude immunized agents
temp = randperm(length(num_list));
rand_infect = num_list(temp(1:ceil(num*init_infe_rate)));
P(rand_infect,8) = latent_period*24*60; % Ip¡úI
P(rand_infect,7) = P(rand_infect,8) - P(rand_infect,6); % S¡úE
P(rand_infect,10) = P(rand_infect,9) + P(rand_infect,8); % I¡úR

% initialize records
ns = zeros(1,6); % number of infected agents in each activity environment at the current moment
ns_t = zeros(simu_weeks*7*24*60,6); % time series of ns
ns_temp = zeros(simu_weeks*7*(end_t - start_t),6); % time series of ns during simulation period everyday
health_status_ratio = zeros(simu_weeks*7*24*60,6); % time series of each health status proportion 
    
%% Simulation     
h = waitbar(0,'Start simulation');
time = start_t; % current time

% loop in each week
for z = 1:simu_weeks
    % loop in each minute
    for i = 1:7*(end_t-start_t)
        time = time + 1;
        
        % internal infection
        [m1, n]=find(P(:,8) ~=0 & P(:,8) > time & P(:,8)-latent_period*24*60 <= time ); % agents of Ip-status
        [m2, n]=find(P(:,8) ~=0 & P(:,8) <= time & P(:,10) ~=0 & P(:,10) > time); % agents of I-status

        p1 = timetable(m1,i); % locations of Ip-status agents
        p2 = timetable(m2,i); % locations of I-status agents
        p = unique([p1;p2]); % remove duplicate locations
        
        % loop in each location
        for j = 1:length(p)
            
            % dining hall
            if ismember(p(j),dining_hall_id) 
                [mm1,nn] = find(timetable(m1,i) == p(j)); % Ip-status agents in dining halls
                [mm2,nn] = find(timetable(m2,i) == p(j)); % I-status agents in dining halls
                mm1 = m1(mm1);
                mm2 = m2(mm2);
                
                [mm, nn] = find(timetable(:,i) == p(j)); % agents in dining_halls
                near_num = min(length(mm)/ capacity_dining_hall(p(j)) * near_num_dining_hall,near_num_dining_hall); % average contacts
                in = setdiff(mm, dining_hall_p{p(j)}) ; % entrance
                out = setdiff(dining_hall_p{p(j)}, mm) ; % exit
                
                % build contact network (in the form of adjacency matrix)
                % remove contacts of agents who have left
                if isempty(dining_hall_near{p(j)}) 
                    dining_hall_near{p(j)} = zeros(length(in));
                else
                    dining_hall_near{p(j)} = [dining_hall_near{p(j)},zeros(size(dining_hall_near{p(j)},1),length(in))];
                    dining_hall_near{p(j)} = [dining_hall_near{p(j)};zeros(length(in),size(dining_hall_near{p(j)},2))];
                end
                dining_hall_p{p(j)} = [dining_hall_p{p(j)};in];
                [temp,temp2,index] = intersect(out,dining_hall_p{p(j)});
                dining_hall_near{p(j)}(index,:) = [];
                dining_hall_near{p(j)}(:,index) = [];
                dining_hall_p{p(j)}(index) = [];
                
                % build new contacts below the upper limit
                if near_num * length(dining_hall_p{p(j)})>=2 && sum(sum(dining_hall_near{p(j)})) < near_num * length(dining_hall_p{p(j)})-2 && isempty(in) == 0
                    
                    % screen out agents of new entrance
                    in_index = length(dining_hall_p{p(j)})-length(in)+1:length(dining_hall_p{p(j)});
                    all_index = 1:length(dining_hall_p{p(j)}); 
                    temp1 = repmat(all_index,1,length(in_index));
                    temp2 = repmat(in_index,length(all_index),1);
                    temp2 = reshape(temp2,size(temp2,1)*size(temp2,2),1);
                    
                    % list possible contacts
                    contact_list = [temp1',temp2];
                    mark1 = (contact_list(:,1) == contact_list(:,2));
                    contact_list(mark1,:)=[];
                    mark2 = (contact_list(:,1) > contact_list(:,2));
                    contact_list(mark2,:)=[contact_list(mark2,2),contact_list(mark2,1)];
                    contact_list = sortrows(contact_list);
                    contact_list = unique(contact_list,'rows');
                    for k = 1:length(all_index)
                        mar = max(ceil(near_num+2 - sum(dining_hall_near{p(j)}(:,k))),0);
                        candidate = find(contact_list(:,1)==k | contact_list(:,2)==k);
                        randp = randperm(length(candidate));
                        contact_list(candidate(randp(1:length(candidate) - mar)),:)=[];
                    end 
                    
                    % build new contacts based on possible contacts and upper limit
                    randlist = randperm(size(contact_list,1));
                    if length(randlist)<= (ceil(near_num * length(dining_hall_p{p(j)})/2) * 2 - sum(sum(dining_hall_near{p(j)})))/2
                        for index = 1:length(randlist)
                            dining_hall_near{p(j)}(contact_list(randlist(index),1),contact_list(randlist(index),2)) = 1;
                            dining_hall_near{p(j)}(contact_list(randlist(index),2),contact_list(randlist(index),1)) = 1;
                        end
                    else
                        for index = 1:(ceil(near_num * length(dining_hall_p{p(j)})/2) * 2 - sum(sum(dining_hall_near{p(j)})))/2
                            dining_hall_near{p(j)}(contact_list(randlist(index),1),contact_list(randlist(index),2)) = 1;
                            dining_hall_near{p(j)}(contact_list(randlist(index),2),contact_list(randlist(index),1)) = 1;
                        end
                    end
                end
                
                % simulation of infection
                for k  = 1:length(dining_hall_p{p(j)})
                    % S-status agents
                    if P(dining_hall_p{p(j)}(k),7)==0 && P(dining_hall_p{p(j)}(k),5)==0
                        p_flag1 = length(intersect(dining_hall_p{p(j)}(dining_hall_near{p(j)}(k,:)==1),mm1));% number of Ip-status neighbors in its contact network
                        p_flag2 = length(intersect(dining_hall_p{p(j)}(dining_hall_near{p(j)}(k,:)==1),mm2));% number of I-status neighbors in its contact network
                        % record infections
                        if rand<(1-(1-p_dining_hall*K)^p_flag1*(1-p_dining_hall)^p_flag2)/P(dining_hall_p{p(j)}(k),11)
                            P(dining_hall_p{p(j)}(k),7) = time;
                            P(dining_hall_p{p(j)}(k),8) = P(dining_hall_p{p(j)}(k),7) + P(dining_hall_p{p(j)}(k),6);
                            P(dining_hall_p{p(j)}(k),10) = P(dining_hall_p{p(j)}(k),9) + P(dining_hall_p{p(j)}(k),8);
                            ns(1)=ns(1)+1;
                        end
                    end
                end
                
            % lecture hall
            elseif ismember(p(j),lecture_hall_id)
                [mm1,nn] = find(timetable(m1,i) == p(j)); % Ip-status agents in lecture halls
                [mm2,nn] = find(timetable(m2,i) == p(j)); % I-status agents in lecture halls
                mm1 = m1(mm1);
                mm2 = m2(mm2);
                
                major_1 = P(mm1,4); % majors of Ip-status agents
                major_2 = P(mm2,4); % majors of I-status agents
                major = unique([major_1;major_2]); % remove duplicate majors
                
                % loop each influenced major
                for k  = 1:length(major)
                    [mm, nn] = find(timetable(:,i) == p(j) & (P(:,4) == major(k))); % all agents of the major
                    near_num = min(length(mm)/ capacity_lecture(major(k)) * near_num_lecture_hall,near_num_lecture_hall);% average contacts
                    in = setdiff(mm, class_p{major(k)}) ; % entrance
                    out = setdiff(class_p{major(k)}, mm) ; % exit
                    
                    % build contact network (in the form of adjacency matrix)
                    % remove contacts of agents who have left
                    if isempty(class_near{major(k)})
                        class_near{major(k)} = zeros(length(in));
                    else
                        class_near{major(k)} = [class_near{major(k)},zeros(size(class_near{major(k)},1),length(in))];
                        class_near{major(k)} = [class_near{major(k)};zeros(length(in),size(class_near{major(k)},2))];
                    end
                    class_p{major(k)} = [class_p{major(k)};in];
                    [temp,temp2,index] = intersect(out,class_p{major(k)});
                    class_near{major(k)}(index,:) = [];
                    class_near{major(k)}(:,index) = [];
                    class_p{major(k)}(index) = [];
                    
                    % build new contacts below the upper limit
                    if near_num * length(class_p{major(k)})>=2 && sum(sum(class_near{major(k)})) < near_num * length(class_p{major(k)})-2 && isempty(in) == 0
                        
                        % screen out agents of new entrance
                        in_index = length(class_p{major(k)})-length(in)+1:length(class_p{major(k)});
                        all_index = 1:length(class_p{major(k)});
                        temp1 = repmat(all_index,1,length(in_index));
                        temp2 = repmat(in_index,length(all_index),1);
                        temp2 = reshape(temp2,size(temp2,1)*size(temp2,2),1);
                        
                        % list possible contacts
                        contact_list = [temp1',temp2];
                        mark1 = (contact_list(:,1) == contact_list(:,2));
                        contact_list(mark1,:)=[];
                        mark2 = (contact_list(:,1) > contact_list(:,2));
                        contact_list(mark2,:)=[contact_list(mark2,2),contact_list(mark2,1)];
                        contact_list = sortrows(contact_list);
                        contact_list = unique(contact_list,'rows');
                        
                        % build new contacts based on possible contacts and upper limit
                        randlist = randperm(size(contact_list,1));
                        index = 0;
                        while sum(sum(class_near{major(k)})) < ceil(near_num * length(class_p{major(k)})/2) * 2 -2 && index < length(randlist)
                            index = index + 1;
                            if sum(class_near{major(k)}(:,contact_list(randlist(index),1)))<near_num+2 && sum(class_near{major(k)}(:,contact_list(randlist(index),2)))<near_num+2
                                class_near{major(k)}(contact_list(randlist(index),1),contact_list(randlist(index),2)) = 1;
                                class_near{major(k)}(contact_list(randlist(index),2),contact_list(randlist(index),1)) = 1;
                            end
                        end
                    end
                    
                    % simulation of infection
                    for m  = 1:length(class_p{major(k)})
                        % S-status agents
                        if P(class_p{major(k)}(m),7)==0 && P(class_p{major(k)}(m),5)==0
                            p_flag1 = length(intersect(class_p{major(k)}(class_near{major(k)}(m,:)==1),mm1));% number of Ip-status neighbors in its contact network
                            p_flag2 = length(intersect(class_p{major(k)}(class_near{major(k)}(m,:)==1),mm2));% number of I-status neighbors in its contact network
                            % record infections
                            if rand<(1-(1-p_lecture*K)^p_flag1*(1-p_lecture)^p_flag2)/P(class_p{major(k)}(m),11)
                                P(class_p{major(k)}(m),7) = time;
                                P(class_p{major(k)}(m),8) = P(class_p{major(k)}(m),7) + P(class_p{major(k)}(m),6);
                                P(class_p{major(k)}(m),10) = P(class_p{major(k)}(m),9) + P(class_p{major(k)}(m),8);
                                ns(2)=ns(2)+1;
                            end
                        end
                    end
                end
                
            % residential building
            elseif ismember(p(j),residential_id)
                [mm1,nn] = find(timetable(m1,i) == p(j)); % Ip-status agents in residential buildings
                [mm2,nn] = find(timetable(m2,i) == p(j)); % Ip-status agents in residential buildings
                mm1 = m1(mm1);
                mm2 = m2(mm2);
                
                residential_1 = P(mm1,2); % rooms of Ip-status agents
                residential_2 = P(mm2,2); % rooms of I-status agents
                residential =  unique([residential_1;residential_2]); % remove duplicate locations
                
                % loop each influenced room (all agents in the same room are contact network neighbors)
                for k  = 1:length(residential)
                    p_flag1 = length(find(P(mm1,2) == residential(k)));% number of Ip-status neighbors in its contact network
                    p_flag2 = length(find(P(mm2,2) == residential(k)));% number of I-status neighbors in its contact network
                    [mm, nn] = find(timetable(:,i) == p(j) & P(:,2) == residential(k)); % all agents of the room
                    
                    % simulation of infection
                    Q = rand(1,length(mm));
                    [m_i, n_i]=find(Q'<(1-(1-p_residential*K)^p_flag1*(1-p_residential)^p_flag2)./P(mm,11));
                    for l = 1:length(m_i)
                        % record infections of S-status agents
                        if P(mm(m_i(l)),7)==0 && P(mm(m_i(l)),5)==0
                            P(mm(m_i(l)),7) = time;
                            P(mm(m_i(l)),8) = P(mm(m_i(l)),7) + P(mm(m_i(l)),6);
                            P(mm(m_i(l)),10) = P(mm(m_i(l)),9) + P(mm(m_i(l)),8);
                            ns(3)=ns(3)+1;
                        end
                    end
                end
                
            % library
            elseif ismember(p(j),library_id)
                [mm1,nn] = find(timetable(m1,i) == p(j)); % Ip-status agents in library
                [mm2,nn] = find(timetable(m2,i) == p(j)); % I-status agents in library
                mm1 = m1(mm1);
                mm2 = m2(mm2);
                
                [mm, nn] = find(timetable(:,i)==p(j)); % agents in library
                in_all = setdiff(mm, [libaray_p{1};libaray_p{2};libaray_p{3};libaray_p{4};libaray_p{5}]) ; % entrance
                out_all = setdiff([libaray_p{1};libaray_p{2};libaray_p{3};libaray_p{4};libaray_p{5}], mm) ;% exit
                rand_floor = randi([1 5],length(in_all),1);
                floor_list = [];
                
                % determine which floor the infected agents on
                for k = 1:libaray_floors             
                    temp = [libaray_p{k};in_all(rand_floor == k)];
                    inter1 = intersect(temp,mm1);
                    inter2 = intersect(temp,mm2);
                    if  isempty(inter2) == 0  || isempty(inter1) == 0
                        floor_list = [floor_list,k];
                    end
                end
                
                % loop in each floor
                for k  = floor_list
                    in = in_all(rand_floor == k) ;
                    [out,temp,index] = intersect(out_all,libaray_p{k});
                    
                    % build contact network (in the form of adjacency matrix)
                    % remove contacts of agents who have left
                    if isempty(libaray_near{k})
                        libaray_near{k} = zeros(length(in));
                    else
                        libaray_near{k} = [libaray_near{k},zeros(size(libaray_near{k},1),length(in))];
                        libaray_near{k} = [libaray_near{k};zeros(length(in),size(libaray_near{k},2))];
                    end
                    libaray_p{k} = [libaray_p{k};in];
                    libaray_near{k}(index,:) = [];
                    libaray_near{k}(:,index) = [];
                    libaray_p{k}(index) = [];
                    
                    % build new contacts below the upper limit
                    near_num = min(length(libaray_p{k})/ libaray_capa_each * near_num_library,near_num_library);
                    if near_num * length(libaray_p{k})>=2 && sum(sum(libaray_near{k})) < near_num * length(libaray_p{k})-2 && isempty(in) == 0 
                        
                        % screen out agents of new entrance
                        in_index = length(libaray_p{k})-length(in)+1:length(libaray_p{k});
                        all_index = 1:length(libaray_p{k}); 
                        temp1 = repmat(all_index,1,length(in_index));
                        temp2 = repmat(in_index,length(all_index),1);
                        temp2 = reshape(temp2,size(temp2,1)*size(temp2,2),1);
                        
                        % list possible contacts
                        contact_list = [temp1',temp2]; 
                        mark1 = (contact_list(:,1) == contact_list(:,2));
                        contact_list(mark1,:)=[];
                        mark2 = (contact_list(:,1) > contact_list(:,2));
                        contact_list(mark2,:)=[contact_list(mark2,2),contact_list(mark2,1)];
                        contact_list = sortrows(contact_list);
                        contact_list = unique(contact_list,'rows');
                        
                        % build new contacts based on possible contacts and upper limit
                        randlist = randperm(size(contact_list,1));
                        index = 0;
                        while sum(sum(libaray_near{k})) < ceil(near_num * length(libaray_p{k})/2) * 2 -2 && index < length(randlist)
                            index = index + 1;
                            if sum(libaray_near{k}(:,contact_list(randlist(index),1)))< near_num+2 && sum(libaray_near{k}(:,contact_list(randlist(index),2)))< near_num+2
                                libaray_near{k}(contact_list(randlist(index),1),contact_list(randlist(index),2)) = 1;
                                libaray_near{k}(contact_list(randlist(index),2),contact_list(randlist(index),1)) = 1;
                            end
                        end
                    end
                    
                    % simulation of infection
                    for m = 1:length(libaray_p{k})
                        % S-status agents
                        if P(libaray_p{k}(m),7)==0 && P(libaray_p{k}(m),5)==0
                            p_flag1 = length(intersect(libaray_p{k}(libaray_near{k}(m,:)==1),mm1));% number of Ip-status neighbors in its contact network
                            p_flag2 = length(intersect(libaray_p{k}(libaray_near{k}(m,:)==1),mm2));% number of I-status neighbors in its contact network
                            % record infections
                            if rand<(1-(1-p_lecture*K)^p_flag1*(1-p_lecture)^p_flag2)/P(libaray_p{k}(m),11)
                                P(libaray_p{k}(m),7) = time;
                                P(libaray_p{k}(m),8) = P(libaray_p{k}(m),7) + P(libaray_p{k}(m),6);
                                P(libaray_p{k}(m),10) = P(libaray_p{k}(m),9) + P(libaray_p{k}(m),8);
                                ns(4)=ns(4)+1;
                            end
                        end
                    end
                end
                
            % road
            else
                % calculate distance between agents on road
                [mm, nn] = find(timetable(:,i)==p(j));
                coor = squeeze(Timetable(mm,i,:));
                distance = pdist2(coor,coor);
                
                % loop agents on road
                for k = 1:length(mm)
                    % S-status agents
                    if P(mm(k),7)==0 && P(mm(k),5)==0
                        related = find(distance(:,k)<=1); % build contact network 
                        related(related==k)=[];
                        p_flag1 = length(intersect(mm(related),mm1));% number of Ip-status neighbors in its contact network
                        p_flag2 = length(intersect(mm(related),mm2));% number of I-status neighbors in its contact network
                        % record infections
                        if rand<(1-(1-p_road*K)^p_flag1*(1-p_road)^p_flag2)/P(mm(k),11)
                            P(mm(k),7) = time;
                            P(mm(k),8) = P(mm(k),7) + P(mm(k),6);
                            P(mm(k),10) = P(mm(k),9) + P(mm(k),8);
                            ns(5)=ns(5)+1;
                        end
                    end
                end
            end
        end
		
		% external infection
        rand_rate = rand(1,num); 
        infect = find(rand_rate'< p_com./P(:,11));
        for l = 1:length(infect)
            % record infections of S-status agents
            if P(infect(l),7)==0 && P(infect(l),5)==0
                P(infect(l),7) = time;
                P(infect(l),8) = P(infect(l),7) + P(infect(l),6);
                P(infect(l),10) = P(infect(l),9) + P(infect(l),8);
                ns(6)=ns(6)+1;
            end
        end
		
        % record cumulative infection in each activity environment in the time series
        ns_temp(i+(z-1)*(end_t-start_t)*7,:) = ns;
        if mod(i,end_t-start_t)==0
            time = time + 24*60-(end_t-start_t);
        end
        
        % current time (week, day, hour, and minute)
        week = z;
        day = floor((time - (week-1)*7*24*60)/(24*60)) + 1;
        hour = floor(((time - (week-1)*7*24*60) - (day-1)*24*60)/60);
        minute = time - (week-1)*7*24*60 - (day-1)*24*60 - hour*60;
        
		str = [num2str((i+(z-1)*(end_t-start_t)*7)/(simu_weeks*(end_t-start_t)*7)*100),'%',' week',num2str(week),' day',num2str(day),' ',num2str(hour),':',num2str(minute)];
        waitbar((i+(z-1)*(end_t-start_t)*7)/(simu_weeks*(end_t-start_t)*7),h,str)
        if isempty(find(P(:,7)==0))==1
            break
        end
    end
    if isempty(find(P(:,7)==0))==1
        break
    end
end

% fill in the missing data in the time series
if i+(z-1)*(end_t-start_t)*7 <simu_weeks*7*(end_t-start_t)
    ns_temp(i+(z-1)*(end_t-start_t)*7+1:end,:) = repmat(ns,simu_weeks*7*(end_t-start_t)-(i+(z-1)*(end_t-start_t)*7),1);
end
for i = 1:simu_weeks*7
    ns_t(24*60*(i-1)+start_t+1:24*60*(i-1)+end_t,:) = ns_temp((end_t-start_t)*(i-1)+1:(end_t-start_t)*i,:);
    if i<simu_weeks*7
        ns_t(24*60*(i-1)+end_t+1:24*60*i+start_t,:) = repmat(ns_temp((end_t-start_t)*i,:),24*60-(end_t-start_t),1);
    else
        ns_t(24*60*(i-1)+end_t+1:end,:) = repmat(ns_temp((end_t-start_t)*i,:),30,1);
    end
end
close(h)

%% calculate each health status proportion
for time = 1:7*24*60*simu_weeks
    V = length(find(P(:,5)==1));
    E = length(find(P(:,7) ~=0 & P(:,7) <=time & P(:,8)-latent_period*24*60 > time));
    Ip = length(find(P(:,8) ~=0 & P(:,8) > time & P(:,8)-latent_period*24*60 <= time));
    I = length(find(P(:,8) ~=0 & P(:,8) <= time & P(:,10) ~=0 & P(:,10) > time));
    R = length(find(P(:,10) ~=0 & P(:,10) <= time));
    S = num - E -I -R - Ip - V;
    health_status_ratio(time,:) = [S,E,Ip,I,R,V];
end
health_status_ratio = health_status_ratio/num;
active_infectious_ratio = health_status_ratio(:,4);
dlmwrite('Active_infectious_cases.csv',active_infectious_ratio,'precision','%.8f')