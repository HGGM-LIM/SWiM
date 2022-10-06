function [LSreord,  ESreord, FAreord, Ireord, Clabreord, allCVals] = isolines_mod(V,F,S,iso)
% ISOLINES compute a list of isolines for a scalar field S defined on the
% mesh (V,F)
%
% [LS,LD,I] = isolines(V,F,S,iso)
%
% Inputs:
%   V      #V by dim list of vertex positions
%   F      #F by 3 list of triangle indices
%   S      #V list of scalar values defined on V
%   iso    #iso list of iso values
% Outputs:
%   LS     #L by dim list of isoline start positions
%   LD     #L by dim list of isoline destination positions
%   I      #L list of indices into iso revealing corresponding iso value
%   ES     #L by 2 list of the edges containing the isoline start positions
%   ED     #L by 2 list of the edges containing the isoline destination positions
%   FA     #L by dim list of faces containing the isolines
%   FAlab  #L list of labels inside the same isoline ()
%
% alternative: tracing isolines should result in smaller plot data (every
% vertex only appears once
%
% Example:
%   iso = linspace(min(S),max(S),10);
%   [LS,LD] = isolines(V,F,S,iso);
%   colormap(jet(numel(iso)-1));
%   tsurf(F,V,'CData',S,fphong);
%   hold on;
%   plot([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]', ...
%     'Color','k','LineWidth',10);
%   hold off;
%

% make sure iso is a ROW vector
iso = iso(:)';
% number of isolines
niso = numel(iso);

% number of domain positions
n = size(V,1);
% number of dimensions
dim = size(V,2);

% number of faces
m = size(F,1);

% Rename for convenience
S1 = S(F(:,1),:);
S2 = S(F(:,2),:);
S3 = S(F(:,3),:);

% t12(i,j) reveals parameter t where isovalue j falls on the line from
% corner 1 to corner 2 of face i
t12 = bsxfun(@rdivide,bsxfun(@minus,iso,S1),S2-S1);
t23 = bsxfun(@rdivide,bsxfun(@minus,iso,S2),S3-S2);
t31 = bsxfun(@rdivide,bsxfun(@minus,iso,S3),S1-S3);

% replace values outside [0,1] with NaNs
t12( (t12<-eps)|(t12>(1+eps)) ) = NaN;
t23( (t23<-eps)|(t23>(1+eps)) ) = NaN;
t31( (t31<-eps)|(t31>(1+eps)) ) = NaN;

% masks for line "parallel" to 12
l12 = ~isnan(t23) & ~isnan(t31);
l23 = ~isnan(t31) & ~isnan(t12);
l31 = ~isnan(t12) & ~isnan(t23);

% find non-zeros (lines) from t23 to t31
[F12,I12,~] = find(l12);
[F23,I23,~] = find(l23);
[F31,I31,~] = find(l31);
% indices directly into t23 and t31 corresponding to F12 and I12
%ti12 = sub2ind(size(l12),F12,I12);
%ti23 = sub2ind(size(l23),F23,I23);
%ti31 = sub2ind(size(l31),F31,I31);
% faster sub2ind
ti12 = F12+(I12-1)*size(l12,1);
ti23 = F23+(I23-1)*size(l23,1);
ti31 = F31+(I31-1)*size(l31,1);

% compute actual position values
LS = [ ...
    ... % average of vertex positions between 2 and 3
    bsxfun(@times,1-t23(ti12), V(F(F12,2),:)) + ...
    bsxfun(@times,  t23(ti12), V(F(F12,3),:)) ...
    ; ... % average of vertex positions between 2 and 3
    bsxfun(@times,1-t31(ti23), V(F(F23,3),:)) + ...
    bsxfun(@times,  t31(ti23), V(F(F23,1),:)) ...
    ;... % average of vertex positions between 2 and 3
    bsxfun(@times,1-t12(ti31), V(F(F31,1),:)) + ...
    bsxfun(@times,  t12(ti31), V(F(F31,2),:)) ...
    ];
LD = [...
    ... % hcat with average of vertex positions between 3 and 1
    bsxfun(@times,1-t31(ti12), V(F(F12,3),:)) + ...
    bsxfun(@times,  t31(ti12), V(F(F12,1),:)) ...
    ;... % hcat with average of vertex positions between 3 and 1
    bsxfun(@times,1-t12(ti23), V(F(F23,1),:)) + ...
    bsxfun(@times,  t12(ti23), V(F(F23,2),:)) ...
    ;... % hcat with average of vertex positions between 3 and 1
    bsxfun(@times,1-t23(ti31), V(F(F31,2),:)) + ...
    bsxfun(@times,  t23(ti31), V(F(F31,3),:)) ...
    ];

% I is just concatenation of each I12
I = [I12;I23;I31];

% Edges connecting each point
ES = [[F(F12,2) F(F12,3) ];[F(F23,3) F(F23,1)];[F(F31,1) F(F31,2)]];
ED = [[F(F12,3) F(F12,1) ];[F(F23,1) F(F23,2)];[F(F31,2) F(F31,3)]];
FA = [F(F12,:);F(F23,:);F(F31,:)];

% Isolines Labels
sts = nonzeros(unique(I));
Niso = length(sts);
LSreord = [0 0 0];
ESreord = [0 0];
FAreord = [0 0 0];
Ireord = 0;
allCVals = 0;
Clabreord = 0;
for i = 1:Niso
    % Selecting all lines belonging to the same isoline
    ind = find(I == sts(i));
    
    isoFaces = FA(ind,:);
    temp = unique(isoFaces) ;
    [~,ord] = ismember(isoFaces,temp);
    
    Surft.SurfData.vertices = V;
    Surft.SurfData.faces = isoFaces;
    [Surft] = Reorg_Surf(Surft);
    graphSu = surface2graph(Surft);
    
    
    %Dulmage-Mendelsohn decomposition
    [~,p,~,r] = dmperm(graphSu);
    comp_sizes = diff(r);
    num_comps = numel(comp_sizes);
    comps = zeros(1,size(graphSu,1));
    comps(r(1:num_comps)) = ones(1,num_comps);
    comps = cumsum(comps);
    comps(p) = comps;
    
    
    cIndex = comps(ord(:,1))';
    Nclust = max(cIndex);
    cont = 0;
    for j = 1:Nclust
        try
            ind2 = find(cIndex == j);
            ESc = ES(ind(ind2),:); % All the edges (including LS) for the same cluster for the same isoline
            EDc = ED(ind(ind2),:); % All the edges (including LD) for the same cluster for the same isoline
            FAc = FA(ind(ind2),:); % All the faces for the same cluster for the same isoline
            allEdges = [ESc;EDc];  % Concatenating the edges
            allFaces = [FAc;FAc];  % Concatenating the faces (repeating)
            
            [allEdgesS,a,~] = unique(sort(allEdges')','rows'); % Removing duplicate edges (caused by the LS and LD variables)
            Ne = size(allEdgesS,1); % Number of non-repeated edges
            
            allEdgesFAc = [FAc(:,1:2);FAc(:,2:3);FAc(:,[1 3])]; % All the edges from each face
            [~, ordFAc] = ismember(sort(allEdgesFAc(:,1:2)')',allEdgesS,'rows'); % Detecting the edges of the faces which coincide with the edges of the line
            connEdge = sort(reshape(ordFAc,[Ne 3])')';  % Reordering to remove the edge of the face that don not belong to the line
            connEdge = connEdge(:,2:end); % Removing the edge
            
            edgGraph = sparse([connEdge(:,1) connEdge(:,2)],[connEdge(:,2) connEdge(:,1)],1,Ne,Ne);

            
            % Disconnecting one
            start_ed = connEdge(1,1);
            end_ed = connEdge(1,2);
            edgGraph(start_ed, end_ed) = 0;
            edgGraph(end_ed, start_ed) = 0;
            [~, edgesOrder]=graphshortestpath(edgGraph,start_ed,end_ed);
            reordEdges = allEdgesS(edgesOrder,:);
            
            temp = [LS(ind(ind2),:);LD(ind(ind2),:)];
            nVert = temp(a(edgesOrder),:);
            
            % Curvature estimation
            [S,~] =  curve_smooth(nVert,[[1:size(nVert,1)-1]' [2:size(nVert,1)]'],'FixBoundary',true,'MaxIters',20);
            if find(isnan(S))
                S = nVert;
            end
            if size(unique(S,'rows'),1) == size(S,1)
                cont = cont + 1;
                if size(S,1) > 400
                    [curvValues,~] = Internal_Curvature_Estimation(S,20);
                elseif size(S,1) > 10
                    [curvValues,~] = Internal_Curvature_Estimation(S,4);
                else
                    [curvValues,~] = Internal_Curvature_Estimation(S,1);
                end
                curvValues(isnan(curvValues)) = 0;
                allCVals  = [allCVals;mean(curvValues,2)];
                
                
                LSreord = [LSreord;nVert];
                ESreord = [ESreord;reordEdges];
                FAreord = [FAreord;allFaces(a(edgesOrder),:)];
                Ireord = [Ireord; ones(length(edgesOrder),1)*i];
                Clabreord = [Clabreord; ones(length(edgesOrder),1)*cont];
            end
            
        catch ME
            warning(['Isoline in level: ' num2str(i) ' cluster: ' num2str(j) ' failed.']);
            % disp(ME.message);
        end
    end
end
LSreord(1,:) = [];
ESreord(1,:) = [];
Ireord(1,:) = [];
FAreord(1,:) = [];
Clabreord(1,:) = [];
allCVals(1,:) = [];


% % Detecting the neighbors of each isoline
% sts = nonzeros(unique(Ireord));
% Niso = length(sts);
% pairs = [0 0];
% neighMat = zeros(length(Ireord),2);
% l = 0;
% for i = 2:Niso-1
%     ind1 = find(Ireord == sts(i));
%     l = l + length(ind1);
%     ind2 = find(Ireord == sts(i+1));
%     ind3 = find(Ireord == sts(i-1));
%     
%     
%     k1 = dsearchn(LSreord(ind2,:),LSreord(ind1,:));
%     k2 = dsearchn(LSreord(ind3,:),LSreord(ind1,:));
%     
%     
%     pairs = [pairs;ind1 ind2(k1);ind1 ind3(k2)];
%     
% end
% pairs(1,:) = [];
% 
% fact = 1;
% d = sqrt(sum((LSreord(pairs(:,1),:) - LSreord(pairs(:,2),:)).^2,2));
% 
% [m,s] = normfit(d);
% dth = m + fact*s;
% ind2del = find(d > dth);
% 
% pairs(ind2del,:) = [];
% d(ind2del) = [];
% 
% neighMat(pairs(:,1),1) = pairs(:,2);
% neighMat(pairs(:,2),2) = pairs(:,1);




return;

function [kappa,alpha] = Internal_Curvature_Estimation(xPoints,NeighL);


%% ========================== Main Program ============================= %%


Nv = length(xPoints); % Graph dimension
% dGraph = sparse(Nv,Nv); % Creating an empty sparse matrix dimension

pointOrd = [1:Nv]';

% % First order neighbors (Initial Graph)
% Edges = [[1:Nv-1]' [2:Nv]']; % Edges
% if closbool
%     Edges = [Edges;[Nv 1]];
% end


% Loop around Neighbor levels
neighMat = zeros(Nv,2*NeighL);
kappa = zeros(Nv,NeighL);
alpha = zeros(Nv,NeighL);
for i = 1:NeighL
    neighMat(:,[2*i-1 2*i]) = [pointOrd+i pointOrd-i];
    
    % Adjusting the Neighbors Matrix
    ind = find(neighMat <=0);
    neighMat(ind) = neighMat(ind) + Nv;
    ind = find(neighMat > Nv);
    neighMat(ind) = neighMat(ind) - Nv;
    
    
    
    %         [R,~,k] = circumcenter(xPoints',xPoints( neighMat(:,2*i-1),:)',xPoints( neighMat(:,2*i),:)');
    
    
    
    A = xPoints;
    B = xPoints( neighMat(:,2*i-1),:);
    C = xPoints( neighMat(:,2*i),:);
    
    D = cross(B-A,C-A);
    b = normm(A-C);
    c = normm(A-B);
    if nargout == 1
        a = norm(B-C);     % slightly faster if only R is required
        R = a*b*c/2/norm(D);
        return
    end
    E = cross(D,B-A);
    F = cross(D,C-A);
    G = ([b b b].^2.*E-[c c c].^2.*F)./([normm(D) normm(D) normm(D)]).^2/2;
    M = A + G;
    R = normm(G);  % Radius of curvature
    %   if R == 0
    %     k = G;
    %   else
    %     k = G'/R^2;   % Curvature vector
    %   end
    kappa(:,i) = R.^-1;
    %
    %
    %     evI =  xPoints - xPoints( neighMat(:,2*i-1),:);
    %     evO =  xPoints( neighMat(:,2*i),:) - xPoints;
    %
    
    %     % Exterior angles between consecutive edge vectors
    %     alpha(:,i) = atan2( ...
    %         evI(:,1).*evO(:,2) - evI(:,2).*evO(:,1), ...
    %         evI(:,1).*evO(:,1) + evI(:,2).*evO(:,2) );
    %
    %     lI = sqrt(sum(evI.^2,2));
    %     lO = sqrt(sum(evO.^2,2));
    %     kappa(:,i) = alpha(:,i)./ (0.5 * (lI + lO));
    %
    %
    
    
    
    
    
    
    
    %
    %     P1 = xPoints( neighMat(:,2*i-1),:);
    %     P2 = xPoints;
    %     P3 = xPoints( neighMat(:,2*i),:);
    %
    %     n = cross(P1-P2,P3-P2,2); D = -1*dot(n,P1,2);
    %     t = P2-P1; u = P3-P1; v = P3-P2;
    %     t2 = sum((t.^2)')'; u2 = sum((u.^2)')'; w2 = sum((n.^2)')';
    %     c = P1+(repmat(t2.*dot(u,v,2),[1 size(u,2)]).*u-repmat(u2.*dot(t,v,2),[1 size(t,2)]).*t)./(2*repmat(w2,[1 size(P1,2)])+eps);
    %     kappa(:,i) = 1/2*sqrt(t2.*u2.*dot(v,v,2)./w2+eps);
    %
    %
end
return;
