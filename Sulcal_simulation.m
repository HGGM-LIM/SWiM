function Surfo = Sulcal_simulation
close all

% x is an array of x-values.
% mu is the mean
% sig is the standard deviation 
% amp is the (positive or negative)
% vo is the vertical offset from baseline (positive or negative)

% Simulation 1
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
Ngx = 201; % Number of Points of the gaussian
Nlong_points = 500; % Number of Points along the length
Long_Param = 250;

perc = 5;
vo = 0;

x = linspace(-15,15,Ngx);  % Points of the gaussian
l = linspace(0,Long_Param,Nlong_points+1); % Points along the length (10 cm)

% Pairing
pairsIni = [[1:Ngx]' [Ngx:-1:1]'];

allIndexes= repmat(pairsIni(:)',[Nlong_points+1 1]); % Image indexes
allIndexes = allIndexes + Ngx*repmat([0:Nlong_points]',[1 2*Ngx ]); % Indexes in 4D
initP = allIndexes(:,1:Ngx)';
outP  = allIndexes(:,Ngx:-1:1)';
allPairs = [initP(:) outP(:)];
                
                
                
margL = floor((Nlong_points+1)*perc/200); % Margin

%% Depth variation
maxDepth = 20;
percdec = 60; % 100 No effect over the depth

amps = zeros(1,Nlong_points+1); % Depth representation
amps(margL+1:Nlong_points+1-margL) = maxDepth;    % Maximum Depth 2 cm

wfact = 6;
ind2affect = wfact*margL+1:Nlong_points+1-wfact*margL;

tempVar = cos(0:2*pi/(length(ind2affect)-1):2*pi);
y = rescale(tempVar,0,(1 - percdec/100)*maxDepth);
y = y - max(y);
amps(ind2affect) = amps(ind2affect) + y;

amps = smoothts(amps,'e',60); % Smoothing the the depth according to the length 
amps(Nlong_points:-1:Nlong_points/2+1) = amps(1:Nlong_points/2);




% amps2 = smoothts(flipdim(t,2),'e',80); % Smoothing the the depth according to the length 
% amps = (amps + flipdim(amps,2))/2;

%% Width Representation
maxWidth = 4;
percdec = 20; % 100 No effect over the width

sigs = zeros(1,Nlong_points+1); % Width representation
sigs(margL+1:Nlong_points+1-margL) = maxWidth;    % Maximum width 2 cm
wfact = 6;
ind2affect = wfact*margL+1:Nlong_points+1-wfact*margL;

tempVar = cos(0:2*pi/(length(ind2affect)-1):2*pi);
y = rescale(tempVar,0,(1 - percdec/100)*maxWidth);
y = y - max(y);
sigs(ind2affect) = sigs(ind2affect) + y;
sigs = smoothts(sigs,'e',60); % Smoothing the the depth according to the length 
sigs(Nlong_points:-1:Nlong_points/2+1) = sigs(1:Nlong_points/2);
sigs(end) = 0;
% smooth again
sigs = smoothdata(sigs,'rlowess',60);



%% Curvature along the length
% Center of the sulcus
mus = ones(1,Nlong_points+1)*0; %Centro en 10
xos = cos([0:4*pi/(Nlong_points):4*pi]); % deviation representation
xos = xos*0;    % Maximum deviation =  1 cm (Tunning parameter)
% xos = smoothts(xos,'Gaussian','b',5); % Smoothing the the depth according to the length 


allPoints = [0 0 0];
for i =1:length(l)
    mu = mus(i);
    sig = sigs(i);
    amp = amps(i);
    xo = xos(i);
    
    y = gaus(x,mu,sig,amp,vo);
    % Plot gaussian
%     plot3(x+xo, y*0+l(i),-y+min(y), 'b.', 'LineWidth',3)
    allPoints  = [allPoints;[x(:)+xo, y(:)*0+l(i),-y(:)+min(y)]];
    
   % hold on;
    % Add noise
%     yh = y + randn(size(y))*amp*.10;
%      hold on
%     plot(x, yh, 'ro','markersize', 4)
%     grid on
   % title(sprintf('Guassian with \\mu=%.1f \\sigma=%.1f amp=%.1f vo=%.1f', ...
   %     mu, sig, amp, vo));
    
  
    
end
allPoints(1,:) = [];



%%     Creating the simulated sulcal surface
Surfo.SurfData.vertices = allPoints;

% First face indexes
temp = repmat([2:Ngx-1],[2 1]);
fac1 = [1;temp(:);Ngx];

% Second face indexes
fac2 = repmat([Ngx+1:2*Ngx-1],[2 1]);fac2= fac2(:);

% Third face indexes
fac3 = fac1*0;
fac3(1:2:end) = fac1(1:2:end)+1;
fac3(2:2:end) = fac2(2:2:end)+1;

temp = repmat(Ngx*[1:Nlong_points+1-2],[length(fac1) 1]);
Surfo.SurfData.faces = repmat([fac1 fac2 fac3],[Nlong_points 1]) + repmat([zeros(length(fac1),1);temp(:)],[1 3]);

FV2=smoothpatch(Surfo.SurfData,1,2);
Surfo.SurfData.vertices = FV2.vertices;clear FV2;
allPoints = Surfo.SurfData.vertices;
simWidth = sqrt(sum((allPoints(allPairs(:,1),:) - allPoints(allPairs(:,2),:)).^2,2));




Surfo.Is = simWidth;
% col = Val2colors(simWidth);
% ind = find(Surfo.SurfData.vertices(:,3) >= -0.1);
% Surfo.SurfData.FaceVertexCData = col;
% Surfo.SurfData.FaceVertexCData(ind,:) = ones(length(ind),3);

% Plot_Surf(Surfo);
% a = 1;
%%
