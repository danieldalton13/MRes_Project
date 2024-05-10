%GRAND - Ground Structure Analysis and Design Code.
% Tomas Zegard, Glaucio H Paulino - Version 1.0, Dec-2013

%% === MESH GENERATION LOADS/BCS ==========================================
kappa = 1.0; ColTol = 0.999999;
Cutoff = 0.002; Ng = 50; % Plot: Member Cutoff & Number of plot groups

% --- OPTION 2: STRUCTURED-ORTHOGONAL MESH GENERATION ---------------------
[NODE, ELEM, SUPP, LOAD] = StructDomain(60,20,3,1,'bridge');
Lvl = 6; RestrictDomain = []; % No restriction for box domain

%% === GROUND STRUCTURE METHOD ============================================
PlotPolyMesh(NODE,ELEM,SUPP,LOAD) % Plot the base mesh
[BARS] = GenerateGS(NODE,ELEM,Lvl,RestrictDomain,ColTol); % Generate the GS
Nn = size(NODE,1); Ne = length(ELEM); Nb = size(BARS,1);
[BC] = GetSupports(SUPP);                 % Get reaction nodes
[BT,L] = GetMatrixBT(NODE,BARS,BC,Nn,Nb); % Get equilibrium matrix
[F] = GetVectorF(LOAD,BC,Nn);             % Get nodal force vector

fprintf('Mesh: Elements %d, Nodes %d, Bars %d, Level %d\n',Ne,Nn,Nb,Lvl)
BTBT = [BT -BT]; LL = [L; kappa*L]; sizeBTBT = whos('BTBT'); clear BT L
fprintf('Matrix [BT -BT]: %d x %d in %gMB (%gGB full)\n',...
        length(F),length(LL),sizeBTBT.bytes/2^20,16*(2*Nn)*Nb/2^30)

tic, [S, vol, exitflag] = linprog(LL,[],[],BTBT,F,zeros(2*Nb,1));
fprintf('Objective V = %f\nlinprog CPU time = %g s\n',vol,toc);

S = reshape(S,numel(S)/2,2);  % Separate slack variables
A = S(:,1) + kappa*S(:,2);    % Get cross-sectional areas
N = S(:,1) - S(:,2);          % Get member forces

%% === PLOTTING ===========================================================
PlotGroundStructure(NODE,BARS,A,Cutoff,Ng)
PlotBoundary(ELEM,NODE)

%% === StructDomain FUNCTION =============================================
function [NODE, ELEM, SUPP, LOAD] = StructDomain(Nx, Ny, Lx, Ly, ProblemID)
    % Generate structured-orthogonal domains
    [X, Y] = meshgrid(linspace(0, Lx, Nx+1), linspace(0, Ly, Ny+1));
    NODE = [reshape(X, numel(X), 1) reshape(Y, numel(Y), 1)];
    k = 0; ELEM = cell(Nx*Ny, 1);
    for j = 1:Ny, for i = 1:Nx
            k = k + 1;
            n1 = (i-1)*(Ny+1) + j; n2 = i*(Ny+1) + j;
            ELEM{k} = [n1 n2 n2+1 n1+1];
    end, end

    % Define support and load conditions specifically for a bridge
    SUPP = [1 1 1; Nx*(Ny+1)+1 1 1];
    
    % Identifying the top edge nodes
    TopEdgeNodes = find(NODE(:,2) == max(NODE(:,2)));

    % Finding nodes at specific positions along the top edge
    % Adjust positions based on actual node locations
    TopMidSpanNode = TopEdgeNodes(round(length(TopEdgeNodes)/2)); % Mid-span top node
    TopQuarterSpanNode = TopEdgeNodes(round(length(TopEdgeNodes)/4)); % Quarter-span top node
    TopThreeQuarterSpanNode = TopEdgeNodes(round(3*length(TopEdgeNodes)/4)); % Three-quarter-span top node

    % Adding loads at these positions with different intensities
    LOAD = [TopQuarterSpanNode, 0, -2;  % Double load at quarter-span at top
            TopMidSpanNode, 0, -1;  % Load at mid-span at top
            TopThreeQuarterSpanNode, 0, -1.5];  % 1.5 times load at three-quarter-span at top
end
