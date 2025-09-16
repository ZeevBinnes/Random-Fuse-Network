% find_criticality 
% uses binary search, to find the value of critical voltage or critical
% current. 
% under the critical value, the system survives. 
% above the critical value, the system collapses. 

tic

% basic settings

% CDF of critical currents is x^beta.
% beta = 1 means uniform distribution.
beta = 1;

% model: fixed voltage 'V' or fixed current 'I'.
model = 'I';

% number of resistors. 
% l - number of columns. mist be even. 
% w - number of rows.
l = 100;
w = 100;

% number of results to get
num_results = 3;

% how accurate to be in finding I0 or V0 of criticality. 
precision_binary_search = 10^-3;

% syst stores info about the system, and the results.
syst = struct();

for i=1:num_results

    syst(i).model = model;
    syst(i).l = l;
    syst(i).w = w;
    syst(i).beta = beta;   
    syst(i).Ic = rand(w,l).^(1/beta);
    syst(i).s = struct(); % stores the results

    % initialize lower bounds and upper bounds
    % these values are good for most sizes. for extreme sizes, change.
    lb = 0;
    if model == 'I'
        ub = 0.3*l;
    elseif model == 'V'
        ub = 0.45*w;
    end

    step_counter = 0;
    idx_ub = 0;
    idx_lb = 0;

    % binary search finding criticality
    while (ub-lb) >= precision_binary_search
        step_counter = step_counter + 1;
        IV0 = (ub+lb)/2;

        data = fuse_network_attack(syst(i), IV0);

        if data.num_burned(end) % network collapsed
            ub = IV0;
            idx_ub = step_counter;
        else
            lb = IV0;
            idx_lb = step_counter;
        end

        % store all results, also not final values.
        if model == 'I'
            syst(i).s(step_counter).I0 = IV0;
            syst(i).s(step_counter).V0 = data.IV;
        elseif model == 'V'
            syst(i).s(step_counter).V0 = IV0;
            syst(i).s(step_counter).I0 = data.IV;
        end
        syst(i).s(step_counter).num_burned = data.num_burned;
        syst(i).s(step_counter).list = data.list;
    end

    % save the indecies of the results closest to criticality, 
    % ub from above and lb below. 
    syst(i).idx_ub = idx_ub;
    syst(i).idx_lb = idx_lb;   
        
end

toc