% New_V is a vector of vehicles with certain IDs

% 1) This is needed for Priority assignment
% Search for availale vehicles in the fleet
vn = 0;
for v=1:num_v
    if ALL_V(v).status == 0
        vn = vn + 1;
        New_V(vn) = ALL_V(v);
        if size1_exist == 1
            if ismember(ALL_V(v).ID,size1_assigned(:,1)) % size1 assigned could not exist
                pos = find(size1_assigned(:,1) == ALL_V(v).ID);
                % Here I am taking away vehicles form V_ass when they belong to size1 assignments so that to not
                % have they repeat themselves in the assignment
                size1_assigned(pos,:) = [];
                V_ass(pos) = [];
                R_ass(pos) = [];
            end
        end
    end
end

% When ALL_V.status changes from 1 to 0, the corresponding v-r couple
% of size1 assignment must be destroyied

% If I free a vehicle of size 1 and then I run the 'Prepare_for_next_call',
% those vehicles are doubled. I should avoid this.

% If there are available vehicles I add them to the set of vehicles that
% will be used in the main algorithm. Otherwise V will be [].
if vn>0
    new_vehicles = length([New_V(:).ID]);
    V = [New_V]; % In the case there are some available vehicles
else
    new_vehicles = 0;
    V = []; % In the case there are not vehicles available
end

% 2)
R = [];

% Here I check if there are assigned requests from previous calls. If yes, I also
% update the vehicles that will enter the main algorithm.
% (I should actually take into account earlier calls.
if size1_exist == 1 % In the case I have to take into account size1 trips from previous calls
    R = [R_ass]; % Here I for sure have the set of R_ass of the previous calls and R_priority
    V = [V V_ass]; % Here I add to V=New_V or to V=[] the vehicles that are serving size1 trips
end

% Here I look If there are new tasks to take into account
if new_R_exist == 1
    R = [R new_R];
end

if priority_exist == 1
    R = [R R_priority]
end

% N = length([R(:).ID]);
N = length(R);
% num_v = length([V(:).ID]);
num_v = length(V);
if num_v>0
    max_cap = max([V(:).c]);
else
    max_cap = 0;
end

% If in the main I call 'Priority_assignment', then sets R and V will be again modified by it.
% Otherwise we just maintain the R = [R_ass] and V = [V_ass];


% Priority_assignment will take care of update V and R depending on the cases