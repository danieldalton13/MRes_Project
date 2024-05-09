%DDscript using the GRAND approach

%% === DOMAIN DEFINITIONS ===============================================
% Function Definitions for Domain Geometry and Boundary Conditions
function [x] = EG1Domain(Demand,Arg)
    BdBox = [0 8 0 2];
    switch(Demand)
        case('Dist');  x = DistFnc(Arg,BdBox);
        case('BC');    x = BndryCnds(Arg{:},BdBox);
        case('BdBox'); x = BdBox;
        case('PFix');  x = FixedPoints(BdBox);
    end
end

% Compute Distance Functions
function Dist = DistFnc(P,BdBox)
    Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
end

% Specify Boundary Conditions
function [x] = BndryCnds(Node, Element, BdBox)
    % Calculate distances for support conditions
    LeftBottom = sqrt((Node(:,1)-BdBox(1)).^2 + (Node(:,2)-BdBox(3)).^2);
    [~, LeftBottom] = sort(LeftBottom);
    RightBottom = sqrt((Node(:,1)-BdBox(2)).^2 + (Node(:,2)-BdBox(3)).^2);
    [~, RightBottom] = sort(RightBottom);
    Supp = [LeftBottom(1)  0 0;
            RightBottom(1) 0 0];

    % Calculate positions for load conditions
    MidSpan = sqrt((Node(:,1)-sum(BdBox(1:2))/2).^2 + (Node(:,2)-BdBox(4)).^2);
    [~, MidSpan] = sort(MidSpan);
    QuarterSpan = sqrt((Node(:,1)-(BdBox(1)+0.25*(BdBox(2)-BdBox(1)))).^2 + (Node(:,2)-BdBox(4)).^2);
    [~, QuarterSpan] = sort(QuarterSpan);
    ThreeQuarterSpan = sqrt((Node(:,1)-(BdBox(1)+0.75*(BdBox(2)-BdBox(1)))).^2 + (Node(:,2)-BdBox(4)).^2);
    [~, ThreeQuarterSpan] = sort(ThreeQuarterSpan);

    % Apply loads at specific spans
    Load = [
        QuarterSpan(1), 0, -1;      % Load at 1/4 span at top
        MidSpan(1), 0, -1;          % Load at midspan at top
        ThreeQuarterSpan(1), 0, -1  % Load at 3/4 span at top
    ];

    % Return the support and load conditions
    x = {Supp, Load};
end



% Specify Fixed Points
function [PFix] = FixedPoints(BdBox)
    PFix = [sum(BdBox(1:2))/2 BdBox(3)];
end

%% === MESH GENERATION LOADS/BCS ==========================================
kappa = 1.0; ColTol = 0.999999;
Cutoff = 0.002; Ng = 50; % Plot: Member Cutoff & Number of plot groups

% --- OPTION 1: POLYMESHER MESH GENERATION --------------------------------
addpath PolyMesher/
[NODE, ELEM, SUPP, LOAD] = PolyMesher(@EG1Domain, 600, 30);
Lvl = 5; RestrictDomain = [];
rmpath PolyMesher/

%% === GROUND STRUCTURE METHOD ============================================
PlotPolyMesh(NODE, ELEM, SUPP, LOAD) % Plot the base mesh
[BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol); % Generate the GS
Nn = size(NODE,1); Ne = length(ELEM); Nb = size(BARS,1);
[BC] = GetSupports(SUPP);                 % Get reaction nodes
[BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb); % Get equilibrium matrix
[F] = GetVectorF(LOAD, BC, Nn);             % Get nodal force vector

fprintf('Mesh: Elements %d, Nodes %d, Bars %d, Level %d\n', Ne, Nn, Nb, Lvl)
BTBT = [BT -BT]; LL = [L; kappa*L]; sizeBTBT = whos('BTBT'); clear BT L
fprintf('Matrix [BT -BT]: %d x %d in %gMB (%gGB full)\n',...
        length(F),length(LL),sizeBTBT.bytes/2^20,16*(2*Nn)*Nb/2^30)

tic, [S, vol, exitflag] = linprog(LL,[],[],BTBT,F,zeros(2*Nb,1));
fprintf('Objective V = %f\nlinprog CPU time = %g s\n', vol, toc);

S = reshape(S,numel(S)/2,2);  % Separate slack variables
A = S(:,1) + kappa*S(:,2);    % Get cross-sectional areas
N = S(:,1) - S(:,2);          % Get member forces

%% === PLOTTING ===========================================================
PlotGroundStructure(NODE, BARS, A, Cutoff, Ng)
PlotBoundary(ELEM, NODE)
