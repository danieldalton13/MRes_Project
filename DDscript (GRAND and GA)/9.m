%% Ground Structure Optimization using GA and Linear Programming
% This script optimizes a ground structure problem using GA for node displacement 
% and LP for material usage minimization.

%% Parameters
% Optimization parameters
kappa = 1.0; % Weight factor for material usage
ColTol = 0.999999; % Tolerance for collinearity of bars
Cutoff = 0.05; % Minimum cross-sectional area for plotting
Ng = 50; % Number of groups for color coding bars in plots
TolX = 0; % Tolerance for node movement in x
TolY = 0; % Tolerance for node movement in y
RestrictDomain = []; % No restriction for box domain

% Define the nodes where loads are applied and their respective magnitudes
LoadNodes = [3, 6, 9, 12, 15, 18, 21, 24, 27];
LoadMagnitudes = [1, 1, 1, 1, 1, 1, 1, 1, 1];

% Generate the structured orthogonal domain
[NODE, ELEM, SUPP, LOAD, BARS] = StructDomain(8, 2, 3, 1, 'bridge', LoadNodes, LoadMagnitudes);

% Visualize node numbering
visualizeNodes(NODE, SUPP, LOAD);

% Define movable nodes for optimization
movableNodesX = [5, 11, 17, 23];
movableNodesY = [3, 6, 9, 12, 15, 18, 21, 24, 27];

%% Optimization with GA
numVarsX = length(movableNodesX); % Number of variables to optimize (x-coordinates of movable nodes)
numVarsY = length(movableNodesY); % Number of variables to optimize (y-coordinates of movable nodes)
numVars = numVarsX + numVarsY; % Total number of variables

% Combine lower and upper bounds for x and y coordinates
lb = [NODE(movableNodesX, 1) - TolX; NODE(movableNodesY, 2) - TolY]; % Lower bounds for x and y coordinates
ub = [NODE(movableNodesX, 1) + TolX; NODE(movableNodesY, 2) + TolY]; % Upper bounds for x and y coordinates
options = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', false);

% Define the fitness function for the GA
fitnessFcn = @(xy) evaluateStructure(xy, movableNodesX, movableNodesY, NODE, ELEM, SUPP, LOAD, BARS, ColTol, kappa);

% Run the genetic algorithm to find the best node positions
[bestXY, bestVol] = ga(fitnessFcn, numVars, [], [], [], [], lb, ub, [], options);

% Update the NODE coordinates with the best solution found
NODE(movableNodesX, 1) = bestXY(1:numVarsX);
NODE(movableNodesY, 2) = bestXY(numVarsX+1:end);

%% Ground Structure Method
% Plot the base mesh
PlotPolyMesh(NODE, ELEM, SUPP, LOAD);

% Use the initial bars (BARS) directly from the structured domain
Nn = size(NODE, 1); % Number of nodes
Ne = length(ELEM); % Number of elements
Nb = size(BARS, 1); % Number of bars

% Get reaction nodes (supports)
[BC] = GetSupports(SUPP);

% Get equilibrium matrix and member lengths
[BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb);

% Get nodal force vector
[F] = GetVectorF(LOAD, BC, Nn);

fprintf('Mesh: Elements %d, Nodes %d, Bars %d\n', Ne, Nn, Nb);

% Prepare matrix for linear programming
BTBT = [BT -BT]; 
LL = [L; kappa * L]; 
sizeBTBT = whos('BTBT'); 
clear BT L;
fprintf('Matrix [BT -BT]: %d x %d in %gMB (%gGB full)\n', length(F), length(LL), sizeBTBT.bytes/2^20, 16 * (2 * Nn) * Nb / 2^30);

% Solve the linear programming problem to find optimal cross-sectional areas
tic;
[S, vol, exitflag] = linprog(LL, [], [], BTBT, F, zeros(2 * Nb, 1));
fprintf('Objective V = %f\nlinprog CPU time = %g s\n', vol, toc);

% Separate slack variables and compute final areas and member forces
S = reshape(S, numel(S) / 2, 2); % Separate slack variables
A = S(:, 1) + kappa * S(:, 2); % Get cross-sectional areas
N = S(:, 1) - S(:, 2); % Get member forces

%% Plotting Optimized Structure and Angles
% Plot the final ground structure with cross-sectional areas
figure; % Create a new figure window
PlotGroundStructure(NODE, BARS, A, Cutoff, Ng);

% Plot the boundary elements of the mesh
PlotBoundary(ELEM, NODE);

% Calculate and print connection angles at each node in the optimized structure
calculateAndPrintAnglesOptimized(NODE, BARS);

%% Function Definitions

% Function to visualize node numbering
function visualizeNodes(NODE, SUPP, LOAD)
    figure, hold on, axis equal, axis off;
    plot(NODE(:, 1), NODE(:, 2), 'ko'); 
    text(NODE(:, 1) + 0.1, NODE(:, 2), arrayfun(@num2str, (1:size(NODE, 1))', 'UniformOutput', false)); 
    title('Node Numbering');

    % Highlight movable nodes in red for x and y
    movableNodesX = [5, 11, 17, 23];
    movableNodesY = [3, 6, 9, 12, 15, 18, 21, 24, 27];
    plot(NODE(movableNodesX, 1), NODE(movableNodesX, 2), 'ro', 'MarkerSize', 10); 
    plot(NODE(movableNodesY, 1), NODE(movableNodesY, 2), 'mo', 'MarkerSize', 10); 

    % Highlight nodes movable in both x and y directions in yellow
    movableNodesXY = intersect(movableNodesX, movableNodesY);
    plot(NODE(movableNodesXY, 1), NODE(movableNodesXY, 2), 'yo', 'MarkerSize', 10); 

    % Highlight support nodes in blue
    supportNodes = SUPP(:, 1);
    plot(NODE(supportNodes, 1), NODE(supportNodes, 2), 'bo', 'MarkerSize', 10); 

    % Highlight load nodes in green
    loadNodes = LOAD(:, 1);
    plot(NODE(loadNodes, 1), NODE(loadNodes, 2), 'go', 'MarkerSize', 10); 

    legend({'Nodes', 'Movable in X', 'Movable in Y', 'Movable in X and Y', 'Support Nodes', 'Load Nodes'}, ...
           'Location', 'best', 'TextColor', 'black', 'FontSize', 12);
end

% Function to generate the structured domain
function [NODE, ELEM, SUPP, LOAD, BARS] = StructDomain(Nx, Ny, Lx, Ly, ProblemID, LoadNodes, LoadMagnitudes)
    % Generate structured-orthogonal domains with alternating diagonals
    [X, Y] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));
    NODE = [reshape(X, numel(X), 1) reshape(Y, numel(Y), 1)];
    k = 0; ELEM = cell(2*Nx*Ny, 1); % Increase size to accommodate diagonal elements
    for j = 1:Ny
        for i = 1:Nx
            k = k + 1;
            n1 = (i-1)*(Ny+1) + j; 
            n2 = i*(Ny+1) + j;
            n3 = n2 + 1;
            n4 = n1 + 1;
            ELEM{k} = [n1 n2 n3 n4]; % Add rectangular element

            % Add alternating diagonal elements
            k = k + 1;
            if mod(i+j, 2) == 0
                ELEM{k} = [n1 n3]; % Diagonal from bottom-left to top-right
            else
                ELEM{k} = [n2 n4]; % Diagonal from bottom-right to top-left
            end
        end
    end

    % Define support conditions specifically for a bridge
    SUPP = [1 1 1; Nx*(Ny+1)+1 1 1];
    
    % Apply downward loads to the specified nodes with given magnitudes
    LOAD = arrayfun(@(n, m) [n, 0, -m], LoadNodes, LoadMagnitudes, 'UniformOutput', false);
    LOAD = vertcat(LOAD{:}); % Convert cell array to matrix

    % Generate bars (connections) from elements
    BARS = [];
    for i = 1:length(ELEM)
        elem = ELEM{i};
        if length(elem) == 2 % It's a bar
            BARS = [BARS; elem];
        elseif length(elem) == 4 % It's a rectangle
            BARS = [BARS; elem(1) elem(2); elem(2) elem(3); elem(3) elem(4); elem(4) elem(1)];
        end
    end
end

% Function to plot the mesh
function [] = PlotPolyMesh(NODE, ELEM, SUPP, LOAD)
    figure, hold on, axis equal, axis off;
    MaxNVer = max(cellfun(@numel, ELEM)); % Maximum number of vertices in mesh
    PadWNaN = @(E) [E NaN(1, MaxNVer - numel(E))]; % Pad cells with NaN
    ElemMat = cellfun(PadWNaN, ELEM, 'UniformOutput', false);
    ElemMat = vertcat(ElemMat{:}); % Create padded element matrix
    patch('Faces', ElemMat, 'Vertices', NODE, 'FaceColor', 'w'); % Plot mesh elements

    if (nargin == 4 && ~isempty(SUPP) && ~isempty(LOAD))
        plot(NODE(SUPP(:, 1), 1), NODE(SUPP(:, 1), 2), 'b>', 'MarkerSize', 8); % Plot supports
        plot(NODE(LOAD(:, 1), 1), NODE(LOAD(:, 1), 2), 'm^', 'MarkerSize', 8); % Plot loads
    end
    axis tight, drawnow;
end

% Function to get support boundary conditions
function [BC] = GetSupports(SUPP)
    % Return degrees-of-freedom with fixed (prescribed) displacements
    Nf = sum(sum(~isnan(SUPP(:,2:3))));
    BC = zeros(Nf,1); 
    j = 0;
    for i=1:size(SUPP,1)
        if ~isnan(SUPP(i,2)), j = j + 1; BC(j) = 2*SUPP(i) - 1; end
        if ~isnan(SUPP(i,3)), j = j + 1; BC(j) = 2*SUPP(i); end
    end
    if j~=Nf, error('Parsing number mismatch on BCs.'); end
end

% Function to get equilibrium matrix and member lengths
function [BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb)
    % Generate equilibrium matrix BT and get member lengths L
    D = [NODE(BARS(:,2),1)-NODE(BARS(:,1),1) NODE(BARS(:,2),2)-NODE(BARS(:,1),2)];
    L = sqrt(D(:,1).^2+D(:,2).^2);
    D = [D(:,1)./L D(:,2)./L];
    BT = sparse([2*BARS(:,1)-1 2*BARS(:,1) 2*BARS(:,2)-1 2*BARS(:,2)],...
                 repmat((1:Nb)',1,4),[-D D],2*Nn,Nb);
    BT(BC,:) = []; % Apply constraints
end

% Function to get nodal force vector
function [F] = GetVectorF(LOAD, BC, Nn)
    % Return nodal force vector
    Nl = sum(sum(~isnan(LOAD(:,2:3))));
    F = sparse([],[],[],2*Nn,1,Nl);
    for i=1:size(LOAD,1)
        n = LOAD(i,1);
        if ~isnan(LOAD(i,2)), F(2*n-1) = LOAD(i,2); end
        if ~isnan(LOAD(i,3)), F(2*n) = LOAD(i,3); end
    end
    F(BC) = []; % Apply constraints
end

% Function to plot the final ground structure
function [] = PlotGroundStructure(NODE, BARS, A, Cutoff, Ng)
    figure, hold on, axis equal, axis off; % Removed 'Visible' parameter

    % Define a colormap for plotting
    cmap = colormap(jet(Ng));  % Create a jet colormap with Ng colors
    colorIndices = ceil(Ng * (A / max(A))); % Normalize areas to [0, Ng]
    colorIndices(colorIndices == 0) = 1; % Ensure indices start at 1
    
    A = A / max(A); % Normalize areas to [0,1]
    ind = find(A > Cutoff); % Find indices of bars with area above cutoff
    MyGroup = ceil(Ng * A(ind)); % Group bars by their areas
    Groups = cell(Ng, 1); % Initialize groups
    for i = 1:Ng
        Groups{i} = ind(find(MyGroup == i)); 
    end
    for i = Ng:-1:1 % Plot each group of bars
        if ~isempty(Groups{i})
            XY = [NODE(BARS(Groups{i}, 1), :) NODE(BARS(Groups{i}, 2), :)];
            GroupArea = mean(A(Groups{i})); % Mean area for this group
            % Plot bars with line width proportional to their area
            plot(XY(:, [1 3])', XY(:, [2 4])', 'LineWidth', 5 * sqrt(GroupArea), 'Color', cmap(i, :))
        end
    end
    fprintf('-PLOT- Cutoff %g, Groups %g, Bars plotted %g\n', Cutoff, Ng, length(ind));

    % Add colorbar
    colorbar('Ticks', linspace(0, 1, Ng), 'TickLabels', num2cell(linspace(min(A), max(A), Ng)));
    colormap(jet(Ng));
    caxis([min(A) max(A)]);
    title(colorbar, 'Cross-sectional area');
end

% Function to plot the boundary of the mesh
function [] = PlotBoundary(ELEM, NODE)
    % Get number of nodes, elements, and edges per element
    Nn = size(NODE,1); % Number of nodes
    Ne = length(ELEM); % Number of elements
    NpE = cellfun(@numel,ELEM); % Nodes per element

    % Initialize face connectivity matrix
    FACE = sparse([],[],[],Nn,Nn,sum(NpE));
    for i=1:Ne
        MyFACE = [ELEM{i}; ELEM{i}(2:end) ELEM{i}(1)];
        for j=1:NpE(i)
            if FACE(MyFACE(1,j),MyFACE(2,j))==0 % New edge - Flag it
                FACE(MyFACE(1,j),MyFACE(2,j)) = i;
                FACE(MyFACE(2,j),MyFACE(1,j)) =-i;
            elseif isnan(FACE(MyFACE(1,j),MyFACE(2,j)))
                error(sprintf('Edge [%d %d] found in >2 elements',MyFACE(:,j)))
            else % Edge belongs to 2 elements: inside domain. Lock it.
                FACE(MyFACE(1,j),MyFACE(2,j)) = NaN;
                FACE(MyFACE(2,j),MyFACE(1,j)) = NaN;
            end
        end
    end
    [BOUND(:,1),BOUND(:,2)] = find(FACE>0);
    BOUND(:,3) = FACE(sub2ind(size(FACE),BOUND(:,1),BOUND(:,2)));
    % Plot boundary edges
    plot([NODE(BOUND(:,1),1) NODE(BOUND(:,2),1)]',[NODE(BOUND(:,1),2) NODE(BOUND(:,2),2)]','k')
end

% Function to evaluate the structure for a given node configuration
function vol = evaluateStructure(xy, MovableNodeIndicesX, MovableNodeIndicesY, NODE, ELEM, SUPP, LOAD, BARS, ColTol, kappa)
    % Update node positions based on the genetic algorithm variables
    NODE(MovableNodeIndicesX, 1) = xy(1:length(MovableNodeIndicesX));
    NODE(MovableNodeIndicesY, 2) = xy(length(MovableNodeIndicesX)+1:end);

    % Use the initial bars directly
    % Calculate the equilibrium matrix and force vector
    Nn = size(NODE, 1);
    Nb = size(BARS, 1);
    [BC] = GetSupports(SUPP);
    [BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb);
    [F] = GetVectorF(LOAD, BC, Nn);

    % Solve the linear programming problem
    BTBT = [BT -BT];
    LL = [L; kappa * L];
    [S, vol, exitflag] = linprog(LL, [], [], BTBT, F, zeros(2 * Nb, 1));
    
    % Check if the solution is feasible
    if exitflag ~= 1
        vol = inf; % Infeasible solution
    end
end

% Function to calculate and print the connection angles at each node
function calculateAndPrintAnglesOptimized(NODE, BARS)
    for i = 1:size(NODE, 1)
        connectedBars = find(BARS(:,1) == i | BARS(:,2) == i);
        numConnections = length(connectedBars);
        if numConnections > 1
            nodePos = NODE(i, :);
            angles = zeros(1, numConnections);
            
            % Calculate angles of each connected bar with the positive x-axis
            for j = 1:numConnections
                bar = BARS(connectedBars(j), :);
                otherNode = bar(bar ~= i); % The other node in the bar
                otherPos = NODE(otherNode, :);
                vec = otherPos - nodePos;
                angle = atan2d(vec(2), vec(1)); % Angle with the horizontal
                if angle < 0
                    angle = angle + 360; % Ensure the angle is in [0, 360]
                end
                angles(j) = angle;
            end
            
            % Sort angles in clockwise order
            angles = sort(angles, 'ascend');
            
            % Calculate differences between consecutive angles to get angles between bars
            angleDiffs = diff([angles, angles(1) + 360]);
            
            % Print angles
            fprintf('Node %d: %.2f', i, angleDiffs(1));
            for a = 2:length(angleDiffs)
                fprintf(', %.2f', angleDiffs(a));
            end
            fprintf('\n');
        end
    end
end
