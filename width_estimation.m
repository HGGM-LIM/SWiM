function widthMap = width_estimation(Pial,depthMap,varargin)

%Surf = Read_Surface('/Users/alberto/Proyectos/width/103818/surf/lh.white');
%Pial = Read_Surface('/Users/alberto/Proyectos/width/103818/surf/lh.pial');
%load('/Users/alberto/Proyectos/width/103818_test-depth2mm_smooth100pialcurv_20gyralfilter.mat');
% depthMap = morph.left.depthmap;
%depthMap = load_vtk('/Users/alberto/Proyectos/width/103818/travel_depth.vtk');
%depthMap = depthMap.Is;
%curv_pial = read_cfiles('/Users/alberto/Proyectos/width/103818/surf/lh.curv.pial.T2.two');

% DEFAULTS:
stepSize = 0.2;
simulated = 0;
depthThreshold = 1.5;
maxWidth = 25;


% ASSIGNED:
params_to_variables = containers.Map({'Simulated','DepthStep','DepthThreshold','MaxWidth'},...
    {'simulated','stepSize','depthThreshold','maxWidth'});
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous funtion to use assignin to this workspace
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
        error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end



INITIAL_DEPTH = depthThreshold;
MAX_WIDTH = maxWidth;

tStart = tic;

% Surf = Read_Surface('/Users/alberto/Proyectos/width/103818-V2/surf/lh.white');
% Pial = Read_Surface('/Users/alberto/Proyectos/width/103818-V2/surf/lh.pial');
% load('/Users/alberto/Proyectos/width/103818_retest-depth2mm_smooth100pialcurv_20gyralfilter.mat');
% curv_pial = read_cfiles('/Users/alberto/Proyectos/width/103818-V2/surf/lh.curv');
% curv = Surf.Is;
depthLevels = [INITIAL_DEPTH:stepSize:max(depthMap)]';
[LSreord,  ESreord, ~, I, Clabreord, mCval] = isolines_mod(Pial.SurfData.vertices,Pial.SurfData.faces,depthMap,depthLevels);
%[Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2]=patchcurvature(Pial.SurfData);

ESreord = cast(ESreord,'uint32');
I = cast(I,'uint32');
Clabreord = cast(Clabreord,'uint32');

vertexNormals = per_vertex_normals(Pial.SurfData.vertices,Pial.SurfData.faces);

if simulated
    vertexNormals = -vertexNormals;
end

coords = [0 0 0];
nSegments = 0;
coordsLevel = [];
coordsCluster = [];
endpoints = [];
segLabels = [];
lastIndex = 0;
%origPoints = I.*0;
destPoints = I.*0;
distancesEstimated = mCval.*0;

angleWithNormal = mCval.*0;


isoLevels = unique(I);
Ireord = I;


allEndpoints = zeros(length(mCval),1);
for ii = 1:length(isoLevels)
    isoLevel = isoLevels(ii);
    clusts = unique(Clabreord(I == isoLevel));
    
    disp(['IsoLevel ' num2str(isoLevel) ' of ' num2str(length(isoLevels))]);
    
    for jj = 1:length(clusts)
        
        isoClust = clusts(jj);
        
        ind = find(I == isoLevel & Clabreord == isoClust);
        
        %% Selecting each cluster inside the same isoline
        coordP  = LSreord(ind,:);
        % Smoothing. Just for visualization pourpose
        [S,~] =  curve_smooth(coordP,[[1:size(coordP,1)-1]' [2:size(coordP,1)]'],'FixBoundary',true,'MaxIters',10);
        
        % Curvature vector for the selected cluster
        dcurvMap = mCval(ind);
        
        % Convolution with a Gaussian to smooth the curvature values
        if size(S,1) > 40
            w = gausswin(10);
        else
            w = gausswin(3);
        end
        %
        w = w/sum(w);
        filtCurv = conv(dcurvMap, w,'same');
        
%        surfCurv = mean(curv_pial(ESreord(ind,:)),2);
%        lineDiff = movingslope(movingslope(dcurvMap));
        
        % simplify curve using douglas peuker ramer algorithm, then, select
        % the points with an angle between 0 and 0.6pi as high-curvature
        % points
        [ps,ix] = dpsimplify(S,0.75);
        pss = [ps(size(ps,1),:); ps; ps(1,:)];
        p1 = pss(1:end-2,:) - pss(2:end-1,:);
        p2 = pss(3:end,:) - pss(2:end-1,:);
        curvAngle = acos(dot(p1,p2,2)./(normm(p1).*normm(p2)));
        ix = ix((curvAngle < 0.6*pi) & (curvAngle >= 0));
        
        
        
        % Detecting the clusters of points with high curvature values
        Nv = length(filtCurv);
        edgesMat  = [[1:Nv]' [2:Nv 1]'];
        
        dGraph = sparse(Nv,Nv);
        edg2Mat = edgesMat(find(sum(ismember(edgesMat,ix),2) >= 1),:); % edg2Mat = edgesMat(find(sum(ismember(edgesMat,ix),2) == 2),:);
        ind2con = sub2ind([Nv Nv],[edg2Mat(:,1);edg2Mat(:,2)],[edg2Mat(:,2);edg2Mat(:,1)]);
        dGraph(ind2con) = 1;
        
        LabNet = Label_Graph_Components(double(dGraph),0);
        
        % Selecting a single maximum for each curvature cluster
        cPosLabels = max(full(LabNet));
        sts = unique(nonzeros(cPosLabels));
        Nclust = length(sts);
        ix = zeros(Nclust,1);
        for z = 1:Nclust
            indc = find(cPosLabels == sts(z));
            [~,loc] = max(filtCurv(indc));
            ix(z,:) = indc(loc);
        end
        
        
        %
        %        si detectamos menos de 2 puntos hacemos PCA, máximas dimensiones
        %         y tomamos los 2 puntos más alejados
        if length(ix) < 2
            [~, scorepca] = pca(coordP);
            [~, maxind] = max(scorepca(:,1));
            [~, minind] = min(scorepca(:,1));
            ix = [minind; maxind];
            %coords = [coords; linecoords(maxCurvPoints,:)];
        end
        
        %% correct ix by choosing the highest curv value in -5:5 neigh only in curves large enough
        if size(S,1) < 7
            NFILT = 0;
        else
            NFILT = 5;
        end
        x = -NFILT:NFILT;
        ind2check = ix + x;
        correctLowInd = ind2check <= 0;
        correctHighInd = ind2check > length(ind);
        correctHighInd = correctHighInd.*length(ind);
        correctLowInd = correctLowInd.*length(ind);
        ind2check = ind2check + correctLowInd - correctHighInd;
        [~,cind] = max(filtCurv(ind2check),[],2);
        ix = diag(ind2check(:,cind));
        ix = unique(ix);
        %%
        
        % Saving the endpoints localization. All
        allEndpoints(ind(ix)) = 1;
        
        
        coords = [coords; coordP(ix,:)];
        coordsLevel = [coordsLevel; ones(length(ix),1,'uint32').*isoLevel];
        coordsCluster = [coordsCluster; ones(length(ix),1,'uint32').*isoClust];
        
        
        endpoints = [endpoints; ind(ix);];
        
        NcurPoints = length(ind);
        connCurvGraph = [[[1:NcurPoints-1]' [2:NcurPoints]'];[NcurPoints 1]];
        conn2del = find(sum(ismember(connCurvGraph, ix),2));
        connCurvGraph(conn2del,:) = [];
        compos = connected_components(connCurvGraph);
        compos(ix) = 0;
        
        % hay que arreglar esto, pero por ahora, si el nº de compos es
        % menor que ind, poner todo a ceros
        if length(compos) < length(ind)
            compos = zeros(1,length(ind));
        end
        
        segLabels = [segLabels; compos' + nSegments];
        segLabels(nSegments == segLabels) = 0;
        nSegments = max(segLabels) + 1;
    end
    
    disp('Matching points for the isoLevel...');
    
    % emparejar puntos de igual nivel y distinto segmento
    % descartamos el primero ya que es la etiqueta que se usa para separar
    % entre segmentos
    currLevelSegments = unique(segLabels(lastIndex+1:end));
    currLevelLabels = segLabels(lastIndex+1:end);
    indCurrentLevel = find(I == isoLevel);
    points = LSreord(indCurrentLevel,:);
    
    labels = unique(currLevelLabels);
    
    for i = 1:length(labels)
        segment = labels(i);
        
        if segment == 0
            continue
        end
        
        %         indQuery = segment == currLevelLabels;
        indQuery = find(segment == currLevelLabels);
        indOthers = find(segment ~= currLevelLabels);
        
        
        normDir = (vertexNormals(ESreord(indCurrentLevel(indQuery),1),:) + vertexNormals(ESreord(indCurrentLevel(indQuery),2),:))/2;
        normDir = normDir./normm(normDir);
        otherPoints = points(indOthers,:);
        
        Md1 = createns(otherPoints,'Distance','euclidean');
        
        for m = 1:length(indQuery)
            currentPoint = indQuery(m);
            
            if ~currentPoint
                continue
            end
            
            queryPoint = points(currentPoint,:);
            currentNormal = normDir(m,:);
            
            dirToOthers = otherPoints - repmat(queryPoint,size(otherPoints,1),1);
            dirToOthers = dirToOthers./normm(dirToOthers);
            
            angles = dot(dirToOthers,repmat(currentNormal,size(otherPoints,1),1),2);
            validDest = find(angles > 0);
            
            if isempty(validDest)
                destPoints(indCurrentLevel(currentPoint)) = 0;
                distancesEstimated(indCurrentLevel(currentPoint)) = Inf;
            else
                [k, d] = knnsearch(Md1,queryPoint,'Distance','euclidean','K',30);
                
                % check k is in candidatePoints list
                firstOccurrence = find(ismember(k,validDest), 1,'first');
                if isempty(firstOccurrence)
%                     warning(['Ojo ' num2str(i) ' ' num2str(m)]);
                    destPoints(indCurrentLevel(currentPoint)) = 0;
                    distancesEstimated(indCurrentLevel(currentPoint)) = Inf;
                    continue
                end
                
                
                angleWithNormal(indCurrentLevel(currentPoint)) = angles(k(firstOccurrence));
                dist = d(firstOccurrence);
                destPoints(indCurrentLevel(currentPoint)) = indCurrentLevel(indOthers(k(firstOccurrence)));
                distancesEstimated(indCurrentLevel(currentPoint)) = dist;
                
            end
        end
        
        
        %         queryPoints = points(indQuery,:);
        %         otherPoints = points(indOthers,:);
        %         [k, ~] = knnsearch(otherPoints,queryPoints,'Distance','cityblock');
        %         dist = normm(queryPoints - otherPoints(k,:));
        %
        %         destPoints(indCurrentLevel(indQuery)) = indCurrentLevel(indOthers(k));
        %         distancesEstimated(indCurrentLevel(indQuery)) = dist;
        
    end
    lastIndex = length(segLabels);
    
end
coords(1,:) = [];

% remove the pairs where dist > 5mm
% destPoints(distancesEstimated > 5) = 0;

% comprobar líneas que cortan la superficie o que van en dirección
% contraria a la normal
indMatches = find(destPoints);
src = LSreord(indMatches,:);
dest = LSreord(destPoints(indMatches),:) - src;
dest = dest./normm(dest);


clearvars Clabreord FAreord neighMat normDirs



[surfIntersect,hitDist,~] = ray_mesh_intersect(src,dest,Pial.SurfData.vertices,Pial.SurfData.faces);
%[nhits, hitDist] = ray_mesh_intersect_multi(src,dest,Pial.SurfData.vertices,Pial.SurfData.faces);
surfIntersect = find(((distancesEstimated(indMatches) > (hitDist + 0.5)) & surfIntersect > 0 & surfIntersect < Inf) | distancesEstimated(indMatches) > MAX_WIDTH);
%surfIntersect = find(((distancesEstimated(indMatches) >= hitDist) & (nhits > 0 & nhits < Inf) & ~mod(nhits,2)));

pointsRecalculate = cast(indMatches(surfIntersect),'uint32');
posDestRec = zeros(length(pointsRecalculate)*4000,2,'uint32');
posDestDistances = zeros(length(pointsRecalculate)*4000,1,'single');
indexBeginEnd = zeros(length(pointsRecalculate),2,'uint32');
lastInd = 0;
reverseStr = '';
for i = 1:length(pointsRecalculate)
    currIndex = pointsRecalculate(i);
    currLevel = I(currIndex);
    currSegment = segLabels(currIndex);
    segInLevel = unique(segLabels(I == currLevel));
    % exclude current segment from segment list
    segInLevel(currSegment == segInLevel) = [];
    segInLevel(segInLevel == 0) = [];
    
    indOthers = find(ismember(segLabels,segInLevel));
    queryPoint = LSreord(currIndex,:);
    otherPoints = LSreord(indOthers,:);
    
    src = repmat(queryPoint,size(otherPoints,1),1);
    dirToQuery = otherPoints - src;
    dirToQuery = dirToQuery./normm(dirToQuery);
    normDir = (vertexNormals(ESreord(currIndex,1),:) + vertexNormals(ESreord(currIndex,2),:))/2;
    normDir = normDir./normm(normDir);
    angles = dot(dirToQuery,repmat(normDir,size(otherPoints,1),1),2);
    candidates = indOthers(angles > 0);
    
    if lastInd+length(candidates) > length(posDestRec)
        posDestRec(lastInd+length(candidates)+20000,:) = 0;
        posDestDistances(lastInd+length(candidates)+20000) = 0;
    end
    
    % [k,dist] = dsearchn(queryPoint,otherPoints(angles > 0,:));
    dist = normm(repmat(queryPoint,length(candidates),1) - otherPoints(angles > 0,:));
    
    indexBeginEnd(i,:) = [1+lastInd lastInd+length(candidates)];
    posDestRec(indexBeginEnd(i,1):indexBeginEnd(i,2),:) = [repmat(currIndex,length(candidates),1) candidates];
    posDestDistances(indexBeginEnd(i,1):indexBeginEnd(i,2)) = dist;
    
    lastInd = lastInd+length(candidates);
    
    if mod(i,1000) == 0
        percentDone = 100 * i / length(pointsRecalculate);
        msg = sprintf('Checking intersections with surface: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    %     [~,distSurface,~] = ray_mesh_intersect(src,-dirToQuery,Pial.SurfData.vertices,Pial.SurfData.faces);
    %     [k,dist] = dsearchn(queryPoint,otherPoints);
    %     validPoints = find((dist < distSurface) & (angles > 0));
    %     [~,reestimatedPoint] = min(dist(validPoints));
    %     dist(validPoints(reestimatedPoint))
end
posDestRec(lastInd+1:end,:) = [];
posDestDistances(lastInd+1:end) = [];
msg = sprintf('Checking intersections with surface: 100\n\n');
fprintf([reverseStr, msg]);

dirToQuery = LSreord(posDestRec(:,2),:) - LSreord(posDestRec(:,1),:);
dirToQuery = dirToQuery./normm(dirToQuery);

% split in N_CHUNKS to avoid memory problems
N_CHUNKS = 6;
indexIncrement = ceil(size(posDestRec,1)/N_CHUNKS);
lastIndex = 1;
distSurface = posDestDistances.*0-1;
for m = 1:N_CHUNKS
    first = lastIndex;
    last = lastIndex + indexIncrement;
    if last > length(distSurface)
        last = length(distSurface);
    end
    
    [~,distSurface(first:last),~] = ray_mesh_intersect(LSreord(posDestRec(first:last,1),:),dirToQuery(first:last,:),Pial.SurfData.vertices,Pial.SurfData.faces);
    
    lastIndex = last + 1;
end

indValid = posDestDistances <= distSurface;

reverseStr = '';
for i = 1:size(indexBeginEnd,1)
    first = indexBeginEnd(i,1);
    last = indexBeginEnd(i,2);
    
    currentPointRecalculation = posDestRec(first:last,:);
    currentDestDistances = posDestDistances(first:last);
    
    if length(unique(currentPointRecalculation(:,1))) > 1
        disp('Ojo, algo falla');
    end
    
    localValid = indValid(first:last);
    localIndices = currentPointRecalculation(localValid,2);
    localDistances = currentDestDistances(localValid);
    
    [newVal,newPoint] = min(localDistances);
    if isempty(newPoint) || (newVal > MAX_WIDTH)
        destPoints(unique(currentPointRecalculation(:,1))) = 0;
        distancesEstimated(unique(currentPointRecalculation(:,1))) = Inf;
    else
        destPoints(unique(currentPointRecalculation(:,1))) = localIndices(newPoint);
        distancesEstimated(unique(currentPointRecalculation(:,1))) = newVal;
    end
    
    if mod(i,1000) == 0
        percentDone = 100 * i / size(indexBeginEnd,1);
        msg = sprintf('Recalculating intersecting lines: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
msg = sprintf('Recalculating intersecting lines: 100\n');
fprintf([reverseStr, msg]);

clearvars indValid dirToQuery distSurface posDestRec posDestDistances indexBeginEnd localInidces

% widthmap

widthContribs = cell(length(Pial.SurfData.vertices),1);
pairPoints = [(1:length(destPoints))' destPoints];

reverseStr = '';
for i = 1:length(destPoints)
    
    if destPoints(i) ~= 0
        widthContribs{ESreord(i,1)} = [ widthContribs{ESreord(i,1)}; distancesEstimated(i)];
        %widthContribs{ESreord(i,2)} = [ widthContribs{ESreord(i,2)}; distancesEstimated(i)];
        
        widthContribs{ESreord(destPoints(i),1)} = [ widthContribs{ESreord(destPoints(i),1)}; distancesEstimated(i) ];
        %widthContribs{ESreord(destPoints(i),2)} = [ widthContribs{ESreord(destPoints(i),1)}; distancesEstimated(i) ];
    end
    
    if mod(i,1000) == 0
        percentDone = 100 * i / length(destPoints);
        msg = sprintf('Estimating width map: %3.1f', percentDone); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    
end
msg = sprintf('Estimating width map: 100\n'); %Don't forget this semicolon
fprintf([reverseStr, msg]);

widthmapPre = cellfun(@median, widthContribs);

Md1 = createns(Pial.SurfData.vertices);
[k, ~] = knnsearch(Md1,coords);
widthmapPre(k) = eps;

nanMap = isnan(widthmapPre);

disp('Estimating values for unlabeled regions');
[Trip] = Vert_Neibp(double(Pial.SurfData.faces),size(Pial.SurfData.vertices,1),size(Pial.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);

indz = find(temp == 0);
temp(indz) = 1;
indsts = find(nanMap);
widthmapPre(indsts) = 0;
while ~isempty(indsts)
    temp1 = nanMap(temp);
    temp1(indz) = 0;
    temp1 = sum(temp1(indsts,:),2);
    ind2process = find(Trip(indsts,2) - temp1 ~=0);
    
    Num = widthmapPre(temp);
    Num(indz) = 0;
    Num = Num(indsts(ind2process),:);
    Den = sum(logical(Num),2);
    Num = sum(Num,2);
    widthmapPre(indsts(ind2process)) = Num./Den;
    nanMap(indsts(ind2process)) = 0;
    indsts(ind2process) = [];
end



opts.nSmooth = 1;
widthMap = Geodesic_Map_Smoothing(Pial,widthmapPre,opts);

toc(tStart);