%% 2024 OPhO P34 Numerical Analysis %%
% This code aims to solve problem 34 on the 2024 OPhO through numerical
% analysis. The base concept is building a grid of nodes and updating their
% voltages through the checkerboard method. Then, one can sum up individual
% voltage differences across the board and divide that by the individual
% resistances to find the current.



%% Settings

% This parameter defines how many node exist along each side of a small square. Total nodes = (3*detail)^2.
detail = 32; 

% How many times to iterate each checkerboard loop.
iterations = 1000; 

%% Small resistance values
% These are the resistances of the segments between each node and their
% neighbor. It turns out it does not actaully vary with the detail of the
% numerical analysis because the increase from depth is canceled out by the
% increase in width.

r1 = 1/20/10^5/.002;
r2 = 1/5/10^5/.002;
r3 = 1/10/10^5/.002;

%% Initializing references

% The following voltage grid is what we will manipulate
v_grid = zeros(3*detail, 3* detail);

% Initialize one side to 13 mV. It does not actually matter which side, so
% this might not match the problem statement.
v_grid(:,2*detail+1:3*detail) = .013; 


% This mainly marks which conductance square each node is in.
r_grid = zeros(3*detail, 3*detail);
r_grid(:,1:detail) = r1;
r_grid(1:detail,detail+1:2*detail) = r2;
r_grid(detail+1:2*detail,detail+1:2*detail) = r3;
r_grid(2*detail+1:3*detail, detail+1:2*detail) = r1;
r_grid(:,2*detail+1:3*detail) = r2;


%% Iterating through the grid

for i = 1:iterations
    for j = 1:3*detail
        for k = 1:3*detail
            if k ~= 1 && k ~= 3*detail % This prevents the edges attached to the leads from being changed
                if mod(j + k,2) == 0 % This enforces a checkerboard pattern.
                    v_grid(j,k) = calculateV(j,k,v_grid,r_grid);
                end
            end
        end
    end

    for j = 1:3*detail
        for k = 1:3*detail
            if k ~= 1 && k ~= 3*detail % This prevents the edges attached to the leads from being changed
                if mod(j + k,2) == 1 % This enforces a checkerboard pattern.
                    v_grid(j,k) = calculateV(j,k,v_grid,r_grid);
                end
            end
        end
    end

    disp("Iteration " + int2str(i) + " done!")
end

% The following lines calculate the current through two seperate
% interfaces, so verify convergence.
current1 = sum((v_grid(:,detail/2+1)-v_grid(:,detail/2))./r1);
current2 = sum((v_grid(:,5*detail/2+1)-v_grid(:,5*detail/2))./r2);
disp(current1);
disp(current2);

imagesc(v_grid) % Output a pretty graph.


%% Calculating new voltages

function [newV] = calculateV(x,y,v_grid,r_grid)
    store = zeros(2,4);
    i = 1; 

    r0 = r_grid(x,y);

    try
        store(1,i) = v_grid(x-1,y); % Find the voltage
        store(2,i) = (r_grid(x-1,y)+r0)/2; % Average the resistance between two nodes. Usually just r0, but handles nodes on the edges
        i = i+1;
    catch % This handles nodes that are on the edge of the big square. Might throw ArrayOutOfBounds.
    end

    try
        store(1,i) = v_grid(x,y-1);
        store(2,i) = (r_grid(x,y-1)+r0)/2;
        i = i+1;
    catch 
    end

    try
        store(1,i) = v_grid(x+1,y);
        store(2,i) = (r_grid(x+1,y)+r0)/2;
        i = i+1;
    catch 
    end

    try
        store(1,i) = v_grid(x,y+1);
        store(2,i) = (r_grid(x,y+1)+r0)/2;
        i = i+1;
    catch 
    end

    i = i-1;

    store = store(1:2,1:i); % Memory trimming from edge nodes
    

   % Solving Kirchoff's current law
   
    left = 0;
    right = 0;

    for j = 1:i
        left = left + store(1,j)/store(2,j);
        right = right + 1/store(2,j);
    end

    newV = left/right;
end

