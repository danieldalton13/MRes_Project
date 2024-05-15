%% DDscript using GRAND - Ground Structure Analysis and Design Code
% This code uses LP (GRAND) to solve the ground structure problem, and then
% GA for the nodal displacement

%% === MESH GENERATION LOADS/BCS ==========================================
% Parameters for the optimization
kappa = 1.0; % Weight factor for material usage
ColTol = 0.999999; % Tolerance for collinearity of bars
Cutoff = 0.002; % Minimum cross-sectional area for plotting
Ng = 50; % Number of groups for color coding bars in plots
Tol = 0; % Tolerance for node - amount of +/- movement in x

% --- STRUCTURED-ORTHOGONAL MESH GENERATION -------------------------------
% Define the nodes where loads are applied and their respective magnitudes
LoadNodes = [3, 6, 9, 12, 15, 18, 21, 24, 27];
LoadMagnitudes = [0, 0, 3, 0, 3, 0, 3, 0, 0]; 

[NODE, ELEM, SUPP, LOAD] = StructDomain(8, 2, 3, 1, 'bridge', LoadNodes, LoadMagnitudes);
Lvl = 1; % Using Level 1 with diagonals in mesh
RestrictDomain = []; % No restriction for box domain

% === VISUALIZE NODE NUMBERING ============================================
% Create a figure to visualize node numbering
figure, hold on, axis equal, axis off;
plot(NODE(:, 1), NODE(:, 2), 'ko'); % Plot all nodes
text(NODE(:, 1) + 0.1, NODE(:, 2), arrayfun(@num2str, (1:size(NODE, 1))', 'UniformOutput', false)); % Number the nodes
title('Node Numbering');

% === IDENTIFY MOVABLE NODES =============================================
% Highlight specific nodes that are allowed to move
movableNodes = [5, 11, 17, 23];
plot(NODE(movableNodes, 1), NODE(movableNodes, 2), 'ro', 'MarkerSize', 10); % Highlight movable nodes

%% === OPTIMISATION WITH GA ===============================================
% Define the genetic algorithm parameters
numVars = length(movableNodes); % Number of variables to optimize (x-coordinates of movable nodes)
lb = NODE(movableNodes, 1) - Tol; % Lower bounds for x-coordinates
ub = NODE(movableNodes, 1) + Tol; % Upper bounds for x-coordinates
options = optimoptions('ga', 'PopulationSize', 50, 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);

% Define the fitness function for the GA
fitnessFcn = @(x) evaluateStructure(x, movableNodes, NODE, ELEM, SUPP, LOAD, Lvl, RestrictDomain, ColTol, kappa);

% Run the genetic algorithm to find the best node positions
[bestX, bestVol] = ga(fitnessFcn, numVars, [], [], [], [], lb, ub, [], options);

% Update the NODE coordinates with the best solution found
NODE(movableNodes, 1) = bestX;

%% === GROUND STRUCTURE METHOD ============================================
% Plot the base mesh
PlotPolyMesh(NODE, ELEM, SUPP, LOAD);

% Generate the ground structure (potential bar connections)
[BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol);
Nn = size(NODE, 1); % Number of nodes
Ne = length(ELEM); % Number of elements
Nb = size(BARS, 1); % Number of bars

% Get reaction nodes (supports)
[BC] = GetSupports(SUPP);

% Get equilibrium matrix and member lengths
[BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb);

% Get nodal force vector
[F] = GetVectorF(LOAD, BC, Nn);

fprintf('Mesh: Elements %d, Nodes %d, Bars %d, Level %d\n', Ne, Nn, Nb, Lvl);

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

%% === PLOTTING ===========================================================
% Plot the final ground structure with cross-sectional areas
PlotGroundStructure(NODE, BARS, A, Cutoff, Ng);
% Plot the boundary elements of the mesh
PlotBoundary(ELEM, NODE);

%% === FUNCTION DEFINITIONS ===============================================
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

% Function to generate the structured domain
function [NODE, ELEM, SUPP, LOAD] = StructDomain(Nx, Ny, Lx, Ly, ProblemID, LoadNodes, LoadMagnitudes)
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
end

% Function to generate the ground structure (potential bars)
function [BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol)
    if nargin<5, ColTol=0.9999; end
    if (nargin<4 || isempty(RestrictDomain)), RestrictDomain=@(~,~)[];
    elseif nargin<3, error('Not enough input arguments.'), end

    % Get element connectivity matrix
    Nn = max(cellfun(@max, ELEM)); % Number of nodes
    Ne = length(ELEM); % Number of elements
    A1 = sparse(Nn, Nn);
    for i = 1:Ne, A1(ELEM{i}, ELEM{i}) = true; end
    A1 = A1 - speye(Nn, Nn); 
    An = A1;

    % Level 1 connectivity
    [J, I] = find(An); % Find connected nodes
    BARS = [I J];
    D = [NODE(I,1)-NODE(J,1) NODE(I,2)-NODE(J,2)]; % Distance between nodes
    L = sqrt(D(:,1).^2+D(:,2).^2); % Length of bars
    D = [D(:,1)./L D(:,2)./L]; % Normalized direction

    % Only return bars {i,j} with i<j (no duplicate bars)
    A = sparse(BARS(:,1),BARS(:,2),true,Nn,Nn);
    [J,I] = find(tril(A)); 
    BARS = [I J];
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
    figure('Name','DDscript - GRAND','NumberTitle','off')
    hold on, axis equal, axis off;
    
    % Define a colormap for plotting
    cmap = colormap(jet(Ng));  % Create a jet colormap with Ng colors
    
    A = A/max(A); % Normalize areas to [0,1]
    ind = find(A>Cutoff); % Find indices of bars with area above cutoff
    MyGroup = ceil(Ng*A(ind)); % Group bars by their areas
    Groups = cell(Ng,1); % Initialize groups
    for i=1:Ng, Groups{i} = ind(find(MyGroup==i)); end
    for i=Ng:-1:1 % Plot each group of bars
        if ~isempty(Groups{i})
            XY = [NODE(BARS(Groups{i},1),:) NODE(BARS(Groups{i},2),:)];
            GroupArea = mean(A(Groups{i})); % Mean area for this group
            % Plot bars with line width proportional to their area
            plot(XY(:,[1 3])',XY(:,[2 4])','LineWidth',5*sqrt(GroupArea),'Color',cmap(i,:))
        end
    end
    fprintf('-PLOT- Cutoff %g, Groups %g, Bars plotted %g\n',Cutoff,Ng,length(ind));
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
function vol = evaluateStructure(x, MovableNodeIndices, NODE, ELEM, SUPP, LOAD, Lvl, RestrictDomain, ColTol, kappa)
    % Update node positions based on the genetic algorithm variables
    NODE(MovableNodeIndices, 1) = x;

    % Generate the ground structure
    [BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol);

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
