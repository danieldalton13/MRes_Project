% DDscript using GRAND - Ground Structure Analysis and Design Code.

%% === MESH GENERATION LOADS/BCS ==========================================
kappa = 1.0; ColTol = 0.999999;
Cutoff = 0.002; Ng = 50; % Plot: Member Cutoff & Number of plot groups

% --- OPTION 2: STRUCTURED-ORTHOGONAL MESH GENERATION ---------------------
[NODE, ELEM, SUPP, LOAD] = StructDomain(8, 2, 3, 1, 'bridge');
Lvl = 6; RestrictDomain = []; % No restriction for box domain

%% === GROUND STRUCTURE METHOD ============================================
PlotPolyMesh(NODE, ELEM, SUPP, LOAD); % Plot the base mesh
[BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol); % Generate the GS
Nn = size(NODE, 1); Ne = length(ELEM); Nb = size(BARS, 1);
[BC] = GetSupports(SUPP); % Get reaction nodes
[BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb); % Get equilibrium matrix
[F] = GetVectorF(LOAD, BC, Nn); % Get nodal force vector

fprintf('Mesh: Elements %d, Nodes %d, Bars %d, Level %d\n', Ne, Nn, Nb, Lvl);
BTBT = [BT -BT]; LL = [L; kappa * L]; sizeBTBT = whos('BTBT'); clear BT L;
fprintf('Matrix [BT -BT]: %d x %d in %gMB (%gGB full)\n', length(F), length(LL), sizeBTBT.bytes/2^20, 16 * (2 * Nn) * Nb / 2^30);

tic;
[S, vol, exitflag] = linprog(LL, [], [], BTBT, F, zeros(2 * Nb, 1));
fprintf('Objective V = %f\nlinprog CPU time = %g s\n', vol, toc);

S = reshape(S, numel(S) / 2, 2); % Separate slack variables
A = S(:, 1) + kappa * S(:, 2); % Get cross-sectional areas
N = S(:, 1) - S(:, 2); % Get member forces

%% === PLOTTING ===========================================================
PlotGroundStructure(NODE, BARS, A, Cutoff, Ng);
PlotBoundary(ELEM, NODE);

%% === FUNCTION DEFINITIONS ===============================================
function [] = PlotPolyMesh(NODE, ELEM, SUPP, LOAD)
    figure, hold on, axis equal, axis off;
    MaxNVer = max(cellfun(@numel, ELEM)); % Max. num. of vertices in mesh
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

    % Apply a uniform downward load to each node along the top edge
    LOAD = arrayfun(@(n) [n, 0, -1], TopEdgeNodes, 'UniformOutput', false);
    LOAD = vertcat(LOAD{:}); % Convert cell array to matrix
end

function [BARS] = GenerateGS(NODE, ELEM, Lvl, RestrictDomain, ColTol)
    % Function code for GenerateGS...
    if nargin<5, ColTol=0.9999; end
    if (nargin<4 || isempty(RestrictDomain)), RestrictDomain=@(~,~)[];
    elseif nargin<3, error('Not enough input arguments.'), end

    % Get element connectivity matrix
    Nn = max(cellfun(@max, ELEM)); Ne = length(ELEM);
    A1 = sparse(Nn, Nn);
    for i = 1:Ne, A1(ELEM{i}, ELEM{i}) = true; end
    A1 = A1 - speye(Nn, Nn); An = A1;

    % Level 1 connectivity
    [J, I] = find(An); % Reversed because find returns values column-major
    BARS = [I J];
    D = [NODE(I,1)-NODE(J,1) NODE(I,2)-NODE(J,2)];
    L = sqrt(D(:,1).^2+D(:,2).^2);  % Length of bars
    D = [D(:,1)./L D(:,2)./L];      % Normalized direction

    % Levels 2 and above
    for i=2:Lvl
        Aold = An; An = logical(An*A1); Gn = An - Aold; % Get NEW bars @ level 'n'
        [J,I] = find(Gn-diag(diag(Gn)));
        if isempty(J), Lvl = i - 1; fprintf('-INFO- No new bars at Level %g\n',Lvl); break, end

        RemoveFlag = RestrictDomain(NODE,[I J]); % Find and remove bars within restriction zone
        I(RemoveFlag) = []; J(RemoveFlag) = [];

        newD = [NODE(I,1)-NODE(J,1) NODE(I,2)-NODE(J,2)];
        L = sqrt(newD(:,1).^2+newD(:,2).^2);
        newD = [newD(:,1)./L newD(:,2)./L];

        % Collinearity Check
        p = 1; m = 1; RemoveFlag = zeros(size(I)); Nb = size(BARS,1);
        for j=1:Nn
            % Find I(p:q) - NEW bars starting @ node 'j'
            for p=p:length(I), if I(p)>=j, break, end, end
            for q=p:length(I), if I(q)>j, break, end, end
            if I(q)>j, q = q - 1; end

            if I(p)==j
                % Find BARS(m:n) - OLD bars starting @ node 'j'
                for m=1:Nb, if BARS(m,1)>=j, break, end, end
                for n=m:Nb, if BARS(n,1)>j, break, end, end
                if BARS(n,1)>j, n = n - 1; end

                if BARS(n,1)==j
                    % Dot products of old vs. new bars. If ~collinear: mark
                    C = max(D(m:n,:)*newD(p:q,:)',[],1);
                    RemoveFlag(p-1+find(C>ColTol)) = true;
                end
            end
        end

        % Remove collinear bars and make symmetric again. Bars that have one
        % angle marked as collinear but the other not, will be spared
        ind = find(RemoveFlag==false);
        H = sparse(I(ind),J(ind),true,Nn,Nn,length(ind));
        [J,I] = find(H+H');
        fprintf('Lvl %2g - Collinear bars removed: %g\n',i,(length(RemoveFlag)-length(I))/2);

        BARS = sortrows([BARS; I J]);
        D = [NODE(BARS(:,1),1)-NODE(BARS(:,2),1) NODE(BARS(:,1),2)-NODE(BARS(:,2),2)];
        L = sqrt(D(:,1).^2+D(:,2).^2);  % Length of bars
        D = [D(:,1)./L D(:,2)./L];      % Normalized direction
    end

    % Only return bars {i,j} with i<j (no duplicate bars)
    A = sparse(BARS(:,1),BARS(:,2),true,Nn,Nn);
    [J,I] = find(tril(A)); BARS = [I J];
end

function [BC] = GetSupports(SUPP)
    % Return degrees-of-freedom with fixed (prescribed) displacements
    Nf = sum(sum(~isnan(SUPP(:,2:3))));
    BC = zeros(Nf,1); j = 0;
    for i=1:size(SUPP,1)
        if ~isnan(SUPP(i,2)), j = j + 1; BC(j) = 2*SUPP(i) - 1; end
        if ~isnan(SUPP(i,3)), j = j + 1; BC(j) = 2*SUPP(i);     end
    end
    if j~=Nf, error('Parsing number mismatch on BCs.'), end
end

function [BT, L] = GetMatrixBT(NODE, BARS, BC, Nn, Nb)
    % Generate equilibrium matrix BT and get member lengths L
    D = [NODE(BARS(:,2),1)-NODE(BARS(:,1),1) NODE(BARS(:,2),2)-NODE(BARS(:,1),2)];
    L = sqrt(D(:,1).^2+D(:,2).^2);
    D = [D(:,1)./L D(:,2)./L];
    BT = sparse([2*BARS(:,1)-1 2*BARS(:,1) 2*BARS(:,2)-1 2*BARS(:,2)],...
                 repmat((1:Nb)',1,4),[-D D],2*Nn,Nb);
    BT(BC,:) = []; % Apply constraints
end

function [F] = GetVectorF(LOAD, BC, Nn)
    % Return nodal force vector
    Nl = sum(sum(~isnan(LOAD(:,2:3))));
    F = sparse([],[],[],2*Nn,1,Nl);
    for i=1:size(LOAD,1)
        n = LOAD(i,1);
        if ~isnan(LOAD(i,2)), F(2*n-1) = LOAD(i,2); end
        if ~isnan(LOAD(i,3)), F(2*n) = LOAD(i,3);   end
    end
    F(BC) = []; % Apply constraints
end

function [] = PlotGroundStructure(NODE, BARS, A, Cutoff, Ng)
    figure('Name','DDscript - GRAND','NumberTitle','off')
    hold on, axis equal, axis off;
    
    % Define a colormap and store it in variable `cmap`
    cmap = colormap(jet(Ng));  % This creates a jet colormap with Ng colors
    
    A = A/max(A); % Normalize to [0,1] areas
    ind = find(A>Cutoff);
    MyGroup = ceil(Ng*A(ind)); % Round up to the closest group of bars
    Groups = cell(Ng,1);       % Store the indices of similar bars
    for i=1:Ng, Groups{i} = ind(find(MyGroup==i)); end
    for i=Ng:-1:1 % Plot each group of similar bars in a single plot call
        if ~isempty(Groups{i})
            XY = [NODE(BARS(Groups{i},1),:) NODE(BARS(Groups{i},2),:)];
            GroupArea = mean(A(Groups{i})); % Mean area for this group
            % Use the colormap for setting line color
            plot(XY(:,[1 3])',XY(:,[2 4])','LineWidth',5*sqrt(GroupArea),'Color',cmap(i,:))
        end
    end
    fprintf('-PLOT- Cutoff %g, Groups %g, Bars plotted %g\n',Cutoff,Ng,length(ind));
end

function [] = PlotBoundary(ELEM, NODE)
    % Get number of nodes, elements, and edges (nodes) per element
    Nn = size(NODE,1); Ne = length(ELEM); NpE = cellfun(@numel,ELEM);

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
    plot([NODE(BOUND(:,1),1) NODE(BOUND(:,2),1)]',[NODE(BOUND(:,1),2) NODE(BOUND(:,2),2)]','k')
end
