Surf = Sulcal_simulation;
Pial = Surf;
depthMap = abs(Surf.SurfData.vertices(:,3));
% curv_pial = discrete_mean_curvature(Surf.SurfData.vertices,Surf.SurfData.faces);

widthSimulated = Surf.Is;
widthEstimated = width_estimation(Surf,depthMap,'Simulated',1,'DepthThreshold',0.01,'DepthStep',0.1,'MaxWidth',40);

diff = widthSimulated - widthEstimated;

