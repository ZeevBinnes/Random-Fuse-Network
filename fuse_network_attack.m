function data = fuse_network_attack(syst_settings, IV0)
%FUSE_NETWORK_ATTACK attack the fuse network
% 
%   data = FUSE_NETWORK_ATTACK(syst_settings, IV0) simulates an attack on a fuse network, and returns data about the attacks results. 
%   SYST_SETTINGS - a struct with the following fields:
%       syst_settings.model - can be 'I' for fixed current, or 'V' for fixed voltage.
%       syst_settings.l - number of columns of resistors (even int)
%       syst_settings.w - number of rows of resistors (int)
%       syst_settings.Ic - critical current for each resistor (l*w matrix)
%       syst_settings.R - resistance of each resistor (l*w matrix).
%           if not given, treated as ones.
%       syst_settings.isRones - true if all R are ones (bool)
%
%   IV0 - the fixed current / voltage (depends on 'model').
%
%   returns:
%   data - a struct with the following fields:
%       data.IV - a vector of the voltages / currents (the not fixed one) for each step in the simulation.
%       data.num_burned - a vector with the number of resistors burned on each step.
%       data.list - a cell of indices of locations of burned resistors for each step. 

    if ~isfield(syst_settings, 'R')
        syst_settings.R = [];
        syst_settings.isRones = true;
    end

    if syst_settings.model == 'I'
        data = current_attack(syst_settings, IV0);
    elseif syst_settings.model == 'V'
        data = voltage_attack(syst_settings, IV0);
    end

%     data.currents ??
%     data.voltages ??
end


function data = current_attack(syst_settings, I0)
%CURRENT_ATTACK attacks the fuse network with a fixed current.
%
%   data = CURRENT_ATTACK(syst_settings, I0)
%   syst_settings - as in FUSE_NETWORK_ATTACK
%   I0 - fixed current, IV0 from FUSE_NETWORK_ATTACK
%   return: data - as in FUSE_NETWORK_ATTACK

    len = 1000; % estimated max num of steps the simulation will run
    num_burned = zeros(len,1);
    V0 = zeros(len,1);
    burned_resistors_list = cell([len,1]);
    counter = 0;
    
    l = syst_settings.l;
    w = syst_settings.w;
    Ic = syst_settings.Ic;
    % I0 = syst_settings.IV0;
    
    [burned_resistors,graph,active_junctions,R_map] = restart_system(w,l);
    
    % initialize the resistance of the resistors.
    if syst_settings.R
        R = syst_settings.R;
        burned_resistors(R==inf) = true;
    else
        R = ones(w,l);
    end
    
    [G,b] = Kirchhoff_matrix(w,l/2,R,R_map);

    % main loop.
    new_burned = burned_resistors;    
    cont = true;
    while cont
        % let the avalanche run
        [G,b] = Kirchhoff_update(G,b,new_burned,R,R_map,w);
        R(new_burned) = inf;
        [V, active_junctions, graph] = voltages(G, b, 1, active_junctions, graph, new_burned, R_map);
        I = Currents(V,1,R,R_map);
        if ~any(any(I))
            break % no current is flowing: the network is detached
        end
        % current fix according to I0
        [I,V] = fix_currents_and_voltages(I,V,I0);
        [new_burned, burned_resistors] = burn(I,Ic,burned_resistors);
        
        counter = counter+1;
        % burned_resistors_list
        [r,c] = find(new_burned);
        burned_resistors_list{counter} = [r,c];
        % num_burned
        num_burned(counter) = nnz(new_burned);
        V0(counter) = calc_V0(I0,V,R_map,R);
        cont = any(new_burned, 'all');
    end
    V0 = V0(1:counter);
    num_burned = num_burned(1:counter);
    burned_resistors_list = burned_resistors_list(1:counter);
    if counter == 0
        num_burned=inf;
    end
    
    data.IV = V0;
    data.num_burned = num_burned;
    data.list = burned_resistors_list;

end


function data = voltage_attack(syst_settings, V0)
%VOLTAGE_ATTACK attacks the fuse network with a fixed voltage.
%
%   data = VOLTAGE_ATTACK(syst_settings, I0)
%   syst_settings - as in FUSE_NETWORK_ATTACK
%   V0 - fixed voltage, IV0 from FUSE_NETWORK_ATTACK
%   return: data - as in FUSE_NETWORK_ATTACK

    len = 500; % estimated max num of steps the simulation will run
    num_burned = zeros(len,1);
    I0 = zeros(len,1);
    burned_resistors_list = cell([len,1]);
    counter = 0;
    
    l = syst_settings.l;
    w = syst_settings.w;
    Ic = syst_settings.Ic;
    % V0 = syst_settings.IV0;
    
    [burned_resistors,graph,active_junctions,R_map] = restart_system(w,l);
    
    if syst_settings.R
        R = syst_settings.R;
        burned_resistors(R==inf) = true;
    else
        R = ones(w,l);
    end
    
    [G,b] = Kirchhoff_matrix(w,l/2,R,R_map);

    % main loop.
    new_burned = burned_resistors;    
    cont = true;
    while cont
        % let the avalanche run
        [G,b] = Kirchhoff_update(G,b,new_burned,R,R_map,w);
        R(new_burned) = inf;
        [V, active_junctions, graph] = voltages(G, b, V0, active_junctions, graph, new_burned, R_map);
        I = Currents(V,V0,R,R_map);
        if ~any(any(I))
            break
        end
        [new_burned, burned_resistors] = burn(I,Ic,burned_resistors);
        
        counter = counter+1;
        % burned_resistors_list
        [r,c] = find(new_burned);
        burned_resistors_list{counter} = [r,c];
        % num_burned
        num_burned(counter) = nnz(new_burned);
        I0(counter) = sum(I(1,:));
        cont = any(new_burned, 'all');
    end
    I0 = I0(1:counter);
    num_burned = num_burned(1:counter);
    burned_resistors_list = burned_resistors_list(1:counter);
    if counter == 0
        num_burned=inf;
    end
    
    data.IV = I0;
    data.num_burned = num_burned;
    data.list = burned_resistors_list;

end

function [burned_resistors,graph,active_junctions, R_map] = restart_system(w,l)
%RESTART_SYSTEM initializes the network
%   RESTART_SYSTEM(w,l) initializes a network of w*l resistors.
%   returns:
%   burned_resistors - a w*l boolean matrix of burned resistors
%   graph - a graph representing the network
%   active_junctions - junctions trough which a current can flow. they are
%       not disconnected from the top and bottom af the network.
%   R_map - a mapping from resistors to junctions they touch.

    up = true; % the direction of the upper-right diamon. always up.
    L = l/2;
    burned_resistors = false([w,l]);
    [graph, R_map] = build_graph(L,w,up);
    active_junctions = true(L*(w-1),1);
end


function [new_burned, burned_resistors] = burn(I,Ic,burned_resistors)
%BURN burns resistors
%
%   BURN(I,Ic,burned_resistors) burns the resistors with too big currents.
%   I - matrix of the currents. 
%   Ic - matrix of crtical currents.
%   burned_resistors - boolean matrix of already burned resistors.
%   
%   returns:
%   new_burned - boolean matrix of resistors burned in this step.
%   burned_resistors - boolean matrix of all burned resistors.

    new_burned = (abs(I)>=Ic & ~burned_resistors);
    burned_resistors = new_burned | burned_resistors;
end


function [I,V] = fix_currents_and_voltages(I,V,I0)
%FIX_CURRENTS_AND_VOLTAGES fixes the currents and voltages
%
%   fix_currents_and_voltages(I,V,I0) fixes the currents and voltages
%   computed in the case of ficed current. there, V and I are computed as
%   if V0=1, and then corrected using this function.
%
%   I,V - the computed currents and voltages of whole network, which are wrong.
%   I0 - the fixed (constant) current
%
%   returns: [I,V] - correct currents and voltages of whole network.

    I_first_row = sum(I(1,:));
    if I_first_row == 0
        I = zeros(size(I));
    else
        fix_factor = I0/I_first_row;
        I = I * fix_factor;
        V = V * fix_factor;
    end
end


function V0 = calc_V0(I0,V,R_map,R)  
%CALC_V0 calculates the voltage on the system in total
%
%   CALC_V0(I0,V,R_map,R) calculates the voltage on the system in total
%   I0 - the fixed current
%   V - all voltages on junctions in system
%   R_map - mapping between resistors and junctions
%   R - resistance of all resistors
%
%   returns: V0 - the voltage on the system in total

    rmap = R_map(2,:);
    Rs = R(1,:);    
    Vs = V(rmap)';  % only first row of ressitors and junctions
    V0 = (I0 + sum(Vs./Rs)) / sum(1./Rs);
end


function [V, active_junctions, graph] = voltages(G, b, V0, active_junctions, graph, burned, R_map_G)
%VOLTAGES computes voltages of network junctions.
%
%   voltages(G, b, V0, active_junctions, graph, burned, R_map_G)
%
%   G,b - Kirchhoff matrix and vector
%   V0 - voltage on system in total
%   active_junctions - junctions trough which a current can flow. they are
%       not disconnected from the top and bottom af the network.
%   graph - a graph representing the network
%   burned - boolean matrix of all burned resistors
%   R_map_G - extended mapping from resistors to junctions
%
%   returns:
%   V - voltages of all junctions
%   active_junctions, graph - modified  
%
%   to prevent cases of singular matrices, we find connected components in
%   graph, and ignore junctions which are not active, also in the matrix
%   calculation. 

    active_junctions(~full(diag(G))) = false;
    [active_junctions, graph] = find_active_junctions(burned, R_map_G, graph, active_junctions);
    V = zeros(size(b));
    V(active_junctions) = full(G(active_junctions,active_junctions)\(-1*V0*b(active_junctions)));
end


function [active_junctions, graph] = find_active_junctions(burned, R_map_G, graph, active_junctions)
%FIND_ACTIVE_JUNCTIONS findes active junctions
%
%   FIND_ACTIVE_JUNCTIONS(burned, R_map_G, graph, active_junctions) findes
%   active junctions: junctions which are not disconnected to the edges of
%   the system.
%
%   burned, R_map_G, graph, active_junctions - as in VOLTAGES
%
%   returns: active_junctions, graph, as in VOLTAGES
%
%   to find active junctions, we remove the burned resistors from the
%   graph, and then check which components of the graph are connected to
%   the edges of the system. 

    R1 = R_map_G(1:end-1,:);
    R2 = R_map_G(2:end,:);
    s = R1(burned);
    t = R2(burned);
    graph = rmedge(graph, s, t);
    bins = conncomp(graph);
    rr = (bins == bins(end))';  % the last bin is the bin with the edges of the system.
    active_junctions = active_junctions & rr(1:size(active_junctions));
end


function [GG, R_map_G] = build_graph(L,w,up)
%BUILD_GRAPH builds a graph of the network. 
%
%   BUILD_GRAPH(L,w,up) builds a graph of the network with w*l resistors,
%   and in which the top-left diamond is up.
%   in addition, it builds an extended mapping between the resistors and
%   the junctions.
%
%   returns: GG, R_map_G - the graph and the mapping.

    n_junctions = L*(w-1);
    le = n_junctions + 2*L + 1;
    n = 2*L*w + 2*L;
    st = zeros(2, n);
    R_map = Resistors_map(L,w,up);

    R1 = zeros(2,L);
    for m=1:2*L
        st(:,m) = [le;n_junctions+m];
    end
    
    if up
        for j = 1:L
            R1(1,2*j-1) = n_junctions+j;
            R1(1,2*j) = n_junctions+j;
        end
    else
        R1(1,1) = n_junctions+L;
        for j = 1:L-1
            R1(1,2*j) = n_junctions+j;
            R1(1,2*j+1) = n_junctions+j;
        end
        R1(1,2*L) = n_junctions+L;
    end
    if (up && ~mod(w,2)) || (~up && mod(w,2))
        for j = 1:L
            R1(2,2*j-1) = n_junctions+L+j;
            R1(2,2*j) = n_junctions+L+j;
        end
    else
        R1(2,1) = n_junctions+(2*L);
        for j = 1:L-1
            R1(2,2*j) = n_junctions+L+j;
            R1(2,2*j+1) = n_junctions+L+j;
        end
        R1(2,2*L) = n_junctions+(2*L);
    end
    
    R_map_G = [R1(1,:);R_map;R1(2,:)];
    
    m = 2*L+1;
    for i=1:w
        st(:,m:m+2*L-1) = R_map_G(i:i+1,:);
        m=m+2*L;
    end
    
    GG = graph(st(1,:),st(2,:));
end


function [G,b] = Kirchhoff_matrix(w,L,R,R_map_G)
%Kirchhoff_matrix builds the Kirchhoff equations matrix
%
%   Kirchhoff_matrix(w,L,R,R_map_G) builds a matrix for a system with
%   (w*2L) resistors, with volages R, and mapping R_map_G from resistors to
%   junctions.
%
%   return G,b - Kirchhoff sparse matrix and vector
%
%   the matrix becomes very big. so it is compulsory to use sparse matices. 

    R_map = R_map_G(2:end-1,:);
    
    %build the matrix
    N = (w-1)*L;
    isparse = zeros(1,5*N-(4*L));
    jsparse = zeros(1,5*N-(4*L));
    Gsparse = zeros(1,5*N-(4*L));
    
    d = zeros(N,1);
    for i=1:w-1
        for j=1:2*L
            d(R_map(i,j)) = d(R_map(i,j)) - (1/R(i,j)) - (1/R(i+1,j));
        end
    end
    
    for i = 1:N
        isparse(i) = i;
        jsparse(i) = i;
        Gsparse(i) = d(i);
    end
    
    % V = spdiags(d,0,wL,wL);
    
    m = N;
    
    for i=2:w-1
        for j=1:2*L
            g = 1/R(i,j);
            m = m+1;
            isparse(m) = R_map(i-1,j);
            jsparse(m) = R_map(i,j);
            Gsparse(m) = g;
            m = m+1;
            isparse(m) = R_map(i,j);
            jsparse(m) = R_map(i-1,j);
            Gsparse(m) = g;
        end
    end
    
    G = sparse(isparse,jsparse,Gsparse,N,N);
    % G = sparse(isparse(1:m),jsparse(1:m),Gsparse(1:m),N,N);
    
    bb = zeros(1,L);
    for j=1:2*L
        bb(R_map(1,j)) = bb(R_map(1,j))+1/R(1,j);
    end
    b = sparse(1:L,ones(1,L),bb,N,1);

end


function [G,b] = Kirchhoff_matrix_ones(w,L,R_map_G)
%Kirchhoff_matrix_ones builds the Kirchhoff equations matrix, for the case
%where R==1 everywhere.
%   curruntly not in use

    R_map = R_map_G(2:end-1,:);
    
    %build the matrix
    N = (w-1)*L;
    isparse = zeros(1,5*N-(4*L));
    jsparse = zeros(1,5*N-(4*L));
    Gsparse = zeros(1,5*N-(4*L));

    isparse(1:N) = 1:N;
    jsparse(1:N) = 1:N;
    Gsparse(1:N) = -4;
    
    m = N;
    
    for i=2:w-1
        for j=1:2*L
            m = m+1;
            isparse(m) = R_map(i-1,j);
            jsparse(m) = R_map(i,j);
            Gsparse(m) = 1;
            m = m+1;
            isparse(m) = R_map(i,j);
            jsparse(m) = R_map(i-1,j);
            Gsparse(m) = 1;
        end
    end
    
    G = sparse(isparse,jsparse,Gsparse,N,N);
    % G = sparse(isparse(1:m),jsparse(1:m),Gsparse(1:m),N,N);
    
    bb = zeros(1,L);
    for j=1:2*L
        bb(R_map(1,j)) = bb(R_map(1,j))+1;
    end
    b = sparse(1:L,ones(1,L),bb,N,1);

end


function [G,b] = Kirchhoff_update(G,b,new_burned,R,R_map_G,w)
%Kirchhoff_update updates the Kirchhoff equations matrix.
%
%   Kirchhoff_update(G,b,new_burned,R,R_map_G,w) updates the Kirchhoff
%   matrix according to the new burned resistors, and the previous
%   resistance R they had. 
%
%   return G,b - Kirchhoff sparse matrix and vector
%
%   instead of rebuilding the whoule matrix, we only update it with the new
%   burned resistors, a lot less operations.

    R_map = R_map_G(2:end-1,:);
    
    [r,c] = find(new_burned);
    for n = 1:length(r)
        i = r(n); j = c(n);
        if i == 1
            junc = R_map(i,j);
            G(junc,junc) = diff(G(junc,junc) + 1/R(i,j), 0);
            b(junc) = b(junc) - 1/R(i,j);
        elseif i == w
            junc = R_map(end,j);
            G(junc,junc) = diff(G(junc,junc) + 1/R(end,j), 0);
        else
            junc1 = R_map(i,j);
            junc2 = R_map(i-1,j);
            G(junc1,junc1) = diff(G(junc1,junc1) + 1/R(i,j), 0);
            G(junc2,junc2) = diff(G(junc2,junc2) + 1/R(i,j), 0);
            G(junc1,junc2) = 0;
            G(junc2, junc1) = 0;
        end
    
    end

end


function [G,b] = Kirchhoff_update_ones(G,b,new_burned,R_map_G,w)
%Kirchhoff_update_ones updates the Kirchhoff equations matrix, for the case
%where R==1 everywhere.
%   curruntly not in use
    R_map = R_map_G(2:end-1,:);
    
    [r,c] = find(new_burned);
    for n = 1:length(r)
        i = r(n); j = c(n);
        if i == 1
            junc = R_map(i,j);
            G(junc,junc) = G(junc,junc) + 1;
            b(junc) = b(junc) - 1;
        elseif i == w
            junc = R_map(end,j);
            G(junc,junc) = G(junc,junc) + 1;
        else
            junc1 = R_map(i,j);
            junc2 = R_map(i-1,j);
            G(junc1,junc1) = G(junc1,junc1) + 1;
            G(junc2,junc2) = G(junc2,junc2) + 1;
            G(junc1,junc2) = 0;
            G(junc2, junc1) = 0;
        end
    
    end

end


function I = Currents(V,V0,R,R_map_G)
%Currents calculates the currents in all resistors.
%
%   Currents(V,V0,R,R_map_G) calculates the currents in all resistors,
%   using:
%   V - the voltages on all junctions.
%   V0 - the voltage on the system in total.
%   R - resistance of all resistors.
%   R_map_G - extended mapping from resistors to junctions.
%
%   returns I - the currents on all resistors.

    R_map = R_map_G(2:end-1,:);
    
    I = zeros(size(R));
    for j = 1:size(R_map,2)
        I(1,j) = diff(V0,V(R_map(1,j)))/R(1,j);
    end
    for i = 2:size(R_map,1)
        for j = 1:size(R_map,2)
            I(i,j) = diff(V(R_map(i-1,j)),V(R_map(i,j)))/R(i,j);
        end
    end
    for j = 1:size(R_map,2)
        I(end,j) = diff(V(R_map(end,j)),0)/R(end,j);
    end

end


function c = diff(a,b)
%DIFF rounds numerical errors to 0. 
%
%   diff(a,b) returns 0 if a-b is small enough, or a-b otherwise.

    c=a-b;
    if abs(c) < 10^-10
        c = 0;
    end
end


function R_map = Resistors_map(L, w, up)
%Resistors_map creates a mapping between resistors and junctions.
%
%   Resistors_map(L, w, up) creates a mapping between resistors and junctions.
%   it returns a matrix R_map, of the mapping.
%
%   for every cell in R_map, 
%   R(i,j) and R(i+1,j) both touch the junction in R_map(i,j).

    R_map = zeros(w-1,2*L);
    
    for i=1:w-1
        junc = (i-1)*L;
        if up
            R_map(i,1) = L*i;
            for j = 1:L-1
                R_map(i,2*j) = junc+j;
                R_map(i,2*j+1) = junc+j;
            end
            R_map(i,2*L) = junc+L;
        else
            for j = 1:L
                R_map(i,2*j-1) = junc+j;
                R_map(i,2*j) = junc+j;
            end
        end
        up = ~up;
    end

end

% function data = current_attack_ones(syst_settings)
% 
% len = 1000;
% num_burned = zeros(len,1);
% V0 = zeros(len,1);
% burned_resistors_list = cell([len,1]);
% counter = 0;
% 
% l = syst_settings.L;
% w = syst_settings.w;
% l_junctions = l/2;
% 
% if syst_settings.R
%     R = syst_settings.R;
% else
%     R = false(w,l);
% end
% 
% [graph, R_map] = build_graph(L,w,true);
% active_junctions = true(l_junctions*(w-1),1);
% 
% [KM,b] = Kirchhoff_matrix_ones(w,l_junctions,R_map);
% % main loop.
% % new_burned = burned_resistors;    
% cont = true;
% list_new_burned = [];
% while cont
%     % let the avalanche run
%     [KM,b] = Kirchhoff_update_ones(KM,b,list_new_burned,R_map,w);
%     [V, active_junctions, graph] = voltages(KM, b, 1, active_junctions, graph, list_new_burned, R_map);
%     I = Currents(V,1,R,R_map);
%     if ~any(any(I))
%         break
%     end
%     % current fix according to I0
%     [I,V] = fix_currents_and_voltages(I,V,I0);
%     [new_burned, burned_resistors] = burn(I,Ic,burned_resistors);
%     
%     counter = counter+1;
%     % burned_resistors_list
%     [r,c] = find(new_burned);
%     burned_resistors_list{counter} = [r,c];
%     % num_burned
%     num_burned(counter) = nnz(new_burned);
%     V0(counter) = calc_V0(I0,V,R_map,R);
%     cont = any(new_burned, 'all');
% end
% V0 = V0(1:counter);
% num_burned = num_burned(1:counter);
% burned_resistors_list = burned_resistors_list(1:counter);
% if counter == 0
%     num_burned=inf;
% end
% 
% data.IV = V0;
% data.num_burned = num_burned;
% data.list = burned_resistors_list;
% 
% end
