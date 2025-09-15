function data = fuse_network_attack(syst_settings, IV0)
%     syst_settings.model
%     syst_settings.l
%     syst_settings.w
%     syst_settings.Ic
%     syst_settings.R
%     syst_settings.isRones

    if ~isfield(syst_settings, 'R')
        syst_settings.R = [];
        syst_settings.isRones = true;
    end

    if syst_settings.model == 'I'
        data = current_attack(syst_settings, IV0);
    elseif syst_settings.model == 'V'
        data = voltage_attack(syst_settings, IV0);
    end

%     data.IV
%     data.num_burned
%     data.list
%     data.currents ??
%     data.voltages ??
end


function data = current_attack(syst_settings, I0)

len = 1000;
num_burned = zeros(len,1);
V0 = zeros(len,1);
burned_resistors_list = cell([len,1]);
counter = 0;

l = syst_settings.l;
w = syst_settings.w;
Ic = syst_settings.Ic;
% I0 = syst_settings.IV0;

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
    [V, active_junctions, graph] = voltages(G, b, 1, active_junctions, graph, new_burned, R_map);
    I = Currents(V,1,R,R_map);
    if ~any(any(I))
        break
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

len = 500;
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
    up = true;
    L = l/2;
    burned_resistors = false([w,l]);
    [graph, R_map] = build_graph(L,w,up);
    active_junctions = true(L*(w-1),1);
end

function [new_burned, burned_resistors] = burn(I,Ic,burned_resistors)
    new_burned = (abs(I)>=Ic & ~burned_resistors);
    burned_resistors = new_burned | burned_resistors;
end

function [I,V] = fix_currents_and_voltages(I,V,I0)
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
    rmap = R_map(2,:);
    Rs = R(1,:);
    Vs = V(rmap)';
    V0 = (I0 + sum(Vs./Rs)) / sum(1./Rs);
end

function [V, active_junctions, graph] = voltages(G, b, V0, active_junctions, graph, burned, R_map_G)
    active_junctions(~full(diag(G))) = false;
    [active_junctions, graph] = findconncomp(burned, R_map_G, graph, active_junctions);
    V = zeros(size(b));
    V(active_junctions) = full(G(active_junctions,active_junctions)\(-1*V0*b(active_junctions)));
end

function [active_junctions, graph] = findconncomp(burned, R_map_G, graph, active_junctions)
    R1 = R_map_G(1:end-1,:);
    R2 = R_map_G(2:end,:);
    s = R1(burned);
    t = R2(burned);
    graph = rmedge(graph, s, t);
    bins = conncomp(graph);
    rr = (bins == bins(end))';
    active_junctions = active_junctions & rr(1:size(active_junctions));
end

function [GG, R_map_G] = build_graph(L,w,up)
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
    c=a-b;
    if abs(c) < 10^-10
        c = 0;
    end
end

function R_map = Resistors_map(L, w, up)
% for every cell in R_map, 
% r(i,j) and r(i+1,j) both touch the junction in R_map(i,j).

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
