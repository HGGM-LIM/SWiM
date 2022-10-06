function EstimateSulcalWidth(inputPial,outputWidth,varargin)
%ESTIMATESULCALWIDTH Summary of this function goes here
%   Detailed explanation goes here

depthMap = 1;
outputDepth = 0;

stepSize = 0.2;
depthThreshold = 1.5;
maxWidth = 25;

% ASSIGNED:
params_to_variables = containers.Map({'DepthMap','OutputDepth','DepthStep','DepthThreshold','MaxWidth'},...
    {'depthMap','outputDepth','stepSize','depthThreshold','maxWidth'});
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

Pial = Read_Surface(inputPial);

if depthMap == 1
    depthMap = travelDepth(Pial.SurfData.vertices, Pial.SurfData.faces);
else % if depthMap is a string, load the depth
    depthMap = read_cfiles(depthMap);
end

widthMap = width_estimation(Pial,depthMap,'DepthThreshold',str2double(depthThreshold),'DepthStep',str2double(stepSize),'MaxWidth',str2double(maxWidth));

write_curv(outputWidth, widthMap, length(Pial.SurfData.faces));

if outputDepth ~= 0
    write_curv(outputDepth, depthMap, length(Pial.SurfData.faces));
end

end

