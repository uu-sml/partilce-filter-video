function [] = PF ()
% Just call PF (without any arguments) to run the animation
% 
% This is the matlab code behind the movie "Particle Filter Explained
% without Equations", which can be found at http://youtu.be/aUkBa1zMKv4
% Written by Andreas Svensson, October 2013
% Updated by Andreas Svensson, February 2013, fixing a coding error in the
% 'propagation-update' of the weights
% andreas.svensson@it.uu.se
% http://www.it.uu.se/katalog/andsv164
% 
% The code is provided as is, and I take no responsibility for what this
% code may do to you, your computer or someone else.
%
% This code is licensed under a
% Creative Commons Attribution-ShareAlike 3.0 Unported License.
% http://creativecommons.org/licenses/by-sa/3.0/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup and initialization %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting the random seed, so the same example can be run several times
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% Some unceratinty parameters
measurementNoiseStdev = 0.1; speedStdev = 1;

% Speed of the aircraft
speed = 1;
% Set starting position of aircraft
planePosX = -25; planePosY = 4;

% Some parameters for plotting the particles
m = 1000; k = 0.0001;

% Number of particles
N = 200;

% Some variables for plotting
plotVectorSea = -10:0.01:10;
plotVectorMountains = [-40:0.01:-10.01, 10.01:0.01:40];
plotHeight = 5;

% The function describing the ground
ground = @(x) (x>=10).*((1-(x-10)/30).*sin(x-10)+((x-10)/30).*sin(1.5*(x-10))+0.2.*(x-10).*(x<=20)+2*(x>20))+...
    (x<=-10).*((1-(-x-10)/30).*sin(-x-10)+((-x-10)/30).*sin(1.5*(-x-10))+0.2.*(-x-10).*(x>=-20)+2*(x<-20));

% Plot the environment
area(plotVectorMountains,ground(plotVectorMountains),-1,'FaceColor',[0 0.6 0])
set(gca,'XTick',[]); set(gca,'YTick',[]); hold on
area(plotVectorSea,ground(plotVectorSea),-1,'FaceColor',[0 0 0.8]); axis([-40 40 -1 10])
plane = plotPlane(planePosX,planePosY,1);
measurementLine = line([planePosX planePosX],[ground(planePosX),planePosY],'Color',[1 0 0],'LineStyle',':');
pause(1)


%%%%%%%%%%%%%%%%%%%%%%%
%%% Begin filtering %%%
%%%%%%%%%%%%%%%%%%%%%%%

% Generate particles
particles = rand(N,1)*80-40;

% Plot particles
particleHandle = scatter(particles,plotHeight(ones(size(particles))),m*(1/N*ones(N,1)+k),'k','filled');
pause(1)

FirstRun = 1;

% Initialize particle weights
w = 1/N*ones(N,1);

for t = 1:60
    % Generate height measurements (with gaussian measurement noise)
    planeMeasDist = planePosY - ground(planePosX) + randn*measurementNoiseStdev;
    
    % Evaluate measurements (i.e., create weights) using the pdf for the normal distribution
    w = w.*(1/(sqrt(2*pi)*measurementNoiseStdev)*exp(-((planePosY-ground(particles))-planeMeasDist).^2/(2*measurementNoiseStdev^2)));
    
    % Normalize particle weigths
    w = w/sum(w);

    if FirstRun
        
        % Sort out some particles to evaluate them "in public" the first
        % run (as in the movie)
        [~, order] = sort(w,'descend');
        pmax = order(1);
        pmaxi = setdiff(1:N,pmax);
        delete(particleHandle)
        particleHandle = scatter([particles(pmaxi);particles(pmax)],plotHeight(ones(size(particles))),m*([ones(N-1,1)/N;w(pmax)]+k),'k','filled');
        pause(1)
        
        pmax2 = order(2);
        pmaxi2 = setdiff(pmaxi,pmax2);
        delete(particleHandle)
        particleHandle = scatter([particles(pmaxi2);particles(pmax);particles(pmax2)],plotHeight(ones(size(particles))),m*([ones(N-2,1)/N;w(pmax);w(pmax2)]+k),'k','filled');
        pause(1)
        
        % Plot all weighted particles    
        delete(particleHandle)
        particleHandle = scatter(particles,plotHeight(ones(size(particles))),m*(w+k),'k','filled');
        pause(1)
    end

    % Resample the particles
    u = rand(N,1); wc = cumsum(w);
    [~,ind1] = sort([u;wc]); ind=find(ind1<=N)-(0:N-1)';
    particles=particles(ind,:); w=ones(N,1)./N;

    delete(particleHandle);
    particleHandle = scatter(particles,plotHeight(ones(size(particles))),m*(w+k),'k','filled');
    pause(1)

    % Time propagation
    speedNoise = speedStdev*randn(size(particles));
    particles = particles + speed + speedNoise;
    
    % Update weights
    % w = w, since the update in the previous step is done using our motion model, so the
    % information is already contained in that update.
    
    % Move and plot moved aircraft
    planePosX = planePosX + speed;
    delete(plane); delete(measurementLine)
    plane = plotPlane(planePosX,planePosY,1);
    measurementLine = line([planePosX planePosX],[ground(planePosX),planePosY],'Color',[1 0 0],'LineStyle',':');
    
    if FirstRun
        % Plot updated particles
        delete(particleHandle)
        particleHandle = scatter(particles,plotHeight(ones(size(particles))),m*(w+k),'k','filled');
        pause(1)
    end

    FirstRun = 0;
end
end

function [ h ] = plotPlane( xpos,ypos,fignr )
figure(fignr)

X = xpos - 0.6 + [-1,     -0.1,   -0.09,    0.3,  0.7, 0.8, 0.7, 0.3, -0.09,  -0.1, -1];
Y = ypos + [-0.05, -0.05, -0.4, -0.05, -0.05,0 0.05, 0.05, 0.4, 0.05, 0.05];
h = fill(X,Y,'k');
end