%DDscript using the GRAND approach

%% === MESH GENERATION LOADS/BCS ==========================================
kappa = 1.0; ColTol = 0.999999;
Cutoff = 0.002; Ng = 50; % Plot: Member Cutoff & Number of plot groups

% --- OPTION 2: STRUCTURED-ORTHOGONAL MESH GENERATION ---------------------
[NODE, ELEM, SUPP, LOAD] = StructDomain(60,20,3,1,'MBB');
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

tic, [S,vol,exitflag] = linprog(LL,[],[],BTBT,F,zeros(2*Nb,1));
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

if nargin == 4 || isempty(ProblemID), ProblemID = 1; end
switch ProblemID
    case {'Cantilever','cantilever',1}
        SUPP = [(1:Ny+1)' ones(Ny+1, 2)];
        LOAD = [Nx*(Ny+1)+round((Ny+1)/2) 0 -1];
    case {'MBB','Mbb','mbb',2}
        SUPP = [Nx*(Ny+1)+1 NaN 1;
                (1:Ny+1)' ones(Ny+1, 1) nan(Ny+1, 1)];
        LOAD = [Ny+1 0 -0.5];
    case {'Bridge','bridge',3}
        SUPP = [1 1 1;
                Nx*(Ny+1)+1 1 1];
        LOAD = [(Ny+1)*round(Nx/2)+1 0 -1];
    otherwise
        SUPP = []; LOAD = [];
        disp('-INFO- Structured domain generated with no loads/BC')
end
end
