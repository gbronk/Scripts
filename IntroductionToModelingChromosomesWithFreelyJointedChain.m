% Introduction to Modeling Chromosomes with a Freely Jointed Chain 
% Gabriel Bronk, 2014
% gbronk@brandeis.edu

% This MATLAB code simulates a random conformation of the arm of a chromosome 
% in a yeast cell's nucleus. This script is intended for any researcher who 
% wishes to learn how to model a chromosome in a simple manner. It is also
% useful for teaching college (or advanced high school) students how to
% create a simple Monte Carlo simulation.

% Additional information on this model and more complex versions of this model 
% can be found in the paper: 
% "Effect of Chromosome Tethering on Nuclear Organization in Yeast"  
% by Avsaroglu, Bronk, et al. at:
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102474

% The chromosome is modelled as many straight
% segments, with each segment oriented in a random direction. It's like 
% drawing a zig-zagging line (where each segment is the same length), but
% the line is zig-zagging in 3D. See it now: Chromosomes...in 3D! 
% The nucleus is modelled as a sphere, and the nucleolus is modelled 
% as a spherical cap.

% The code calculates the distance between a particular gene and the 
% spindle pole body (I mean the distance in space, not the
% length of the chromosome). The spindle pole body (SPB) is 
% a structure of proteins located 
% on the nuclear membrane. Many different conformations are simulated, so
% you obtain many different distances between the gene and the SPB. These
% many distances are recorded, and the code plots a probability 
% distribution for the distances.

% Hit the run button to see the code in action before reading through the
% code.

% For information about the commands used in this file:
% Type "help CommandName" into the matlab command window. For instance,  
% type "help clear" (without the quotes), % to find out about the "clear" 
% command. Or googling it gives a lot more info.

clear 
rng('default')
rng('shuffle')
tic

%--------------------------------------------------------------------------

% The following section sets all the parameters in the code. If desired,
% you can change the vaues of these parameters.
NumberOfTrials = 5*10^3; 
% This is how many times to try to simulate chromosome conformations.
% Run this for at least 5*10^5 trials to get a decent sample size.

KuhnLength = 0.200; % units of micrometers
% This is the length of a segment of the chromosome.
ChromosomeCompaction = 13; % units of (micrometers^2)/Mbp   (Mbp is 1 million pase pairs)
% In a given length of a chromosome (in micrometers), there may be more DNA
% (in base pairs) or less DNA. The higher this parameter, the lower the
% amount of DNA in a given length of the chromosome. In this simulation, we
% make the simplifying assumption that chromosome compaction is uniform for
% the entire chromosome.
LengthFromCentromereToTelomere = 0.124; % units of Mbp
% Length of the chromosome arm. 
% The length used here is for the left arm of yeast chromosome III with a
% 0.010 Mbp lac operator array inserted into the arm.
LengthFromCentromereToGene = 0.104; % units of Mbp
% Length of chromosome arm from the centromere to the gene of interest.
% The length used here is the location of the lac operator.
R = 0.95; % units in micrometers
% Radius of the nucleus
Nucleolus = 0.2;
% The fraction of the nucleus' volume that is occupied by the nucleolus.
PeripheryThickness = 0.05; % units in micrometers
% Telomeres are attached to the nuclear envelope. As it is impossible to
% simulate a chromosome whose telomere is exactly at the nuclear envelope,
% we make the aproximation that the telomere can be within a small distance
% of the nuclear envelope. The PeripheryThickness is the maximum distance
% from the nuclear envelope that the telomere can be and still be
% considered attached to the nuclear envelope.
MicrotubuleLength = 0.05; % units in micrometers
% This is the length of the microtubule that connects the chromosome's
% centromere to the SPB.
NumberOfBins = 50;
% Number of bins of the probability distribution.

%--------------------------------------------------------------------------

% This section calculates values based on the parameters above.

SPBCoordinates = [0; 0; R];
% We model the nucleus as a sphere centered at the coordinates (0,0,0).
% We choose the SPB to be located at the very top of the nuclear sphere
% (i.e. coordinates (0, 0, R), and as you will see, the nucleolus will be a
% spherical cap at the bottom of the sphere. This is realistic because in
% real yeast cells, the SPB is at the opposite end of the nucleus from the
% nucleolus.

N = round((ChromosomeCompaction*LengthFromCentromereToTelomere)/((KuhnLength)^2));
% Number of kuhn segments for the whole chromosome arm
n = round((ChromosomeCompaction*LengthFromCentromereToGene)/((KuhnLength)^2));
% Number of kuhn segments between the centromere and the gene of interest.

MaximumDistance = 2*R; 
% The maximum distance to include in the probability distribution
dx = MaximumDistance/NumberOfBins;
% The width of each bin of the probability distribution.
ProbabilityDistribution = zeros(1,NumberOfBins);
% This starts the probabiity distribution at all zeros, and we will
% increase the probabilities in the bins as the simulation progresses.

hroots = roots([(pi/3) (-pi*R) 0 (Nucleolus*4/3*pi*R^3)]);
hindex = find(hroots > 0 & hroots < (2*R));
h = hroots(hindex);
% h is the height of the nucleolus (in micrometers). The lines above do the
% following calculation: We represent the nucleolus as a spherical cap,
% and we want to set this spherical cap to a particular volume:
% e.g. 1/5 of the total nuclear volume. What must the height of the spherical
% cap be in order for the spherical cap to have that particulat volume? The
% first of these three lines solves an algebraic equation that I derived,
% which finds the height of the spherical cap. As the equation is a third
% degree polynomial set equal to 0, the equation has 3 solutions, but only
% one of these solutions makes sense (i.e. is positive and is smaller than
% the diameter of the nuclear sphere). The 2nd and 3rd lines find the
% answer that makes sense.

% -------------------------------------------------------------------------

% Now we simulate the chromosome arm many times. Realistic chromosome
% conformations are those that (1) remain entirely within the nuclear
% sphere, and (2) do not protrude into the nucleolus, and (3) have the telomere 
% in close proximity to the nuclear envelope. We only keep data from simulated 
% conformations that meet all three conditions. 
% I will refer to these conditions later on as "the three necessary conditions".

% The number of times that
% we try to simulate the chromosome arm is "NumberOfTrials". However, most
% of these attempts will fail. This is because most random conformations do not 
% have their telomeres close to the nuclear envelope and/or they protrude 
% out of the nuclear sphere and/or penetrate into the nucleolus.  
% Note that in reality, some protrusion into
% the nucleolus is possible, but our simplified model represents the
% nucleolus as impenetrable),

CoordinateMatrix = zeros(3,N+1);
% This initializes CoordinateMatrix (see below for more info).

for trial = 1:NumberOfTrials
    % The "for loop" tells us to simulate the chromosome arm over and over
    % again until we have simulated it as many times as "NumberOfTrials" is 
    % set to.
    
    CurrentCoordinates = [0; 0; R - MicrotubuleLength];
    % These are the coordinates of the centromere (it's just below the SPB). 
    % As the random walk progresses (i.e. we "draw out" one segment after
    % another), the CurrentCoordinates will be updated to be the
    % coordinates of the end of the segment that has been most recently
    % created.
    CoordinateMatrix(:,1) = [0; 0; R - MicrotubuleLength];
    
    CDFC = norm(CurrentCoordinates); 
    % This is the distance from the CurrentCoordinates to the center of 
    % the nuclar sphere. This stands for "current distance from center".
    u = 0;  
    % u keeps track of the number of segments that we have drawn.
    
    
    while (u < N) && (CDFC < R) && (CurrentCoordinates(3,1) > -R + h)
        % The "while loop" tells MATLAB to keep adding segments until
        % it has added N segments. Also, it tells MATLAB to stop adding 
        % segments if the previous segment has protruded out of the nuclear
        % sphere or into the nucleolus.  If a segment does such protrusion,
        % the while loop does not proceed, the data is not recorded 
        % (i.e. we throw away all the segments of this conformation), and
        % we go to the next trial of the "for loop" in order to start a
        % brand new conformation.
        
        RandomNumbers = randn(3,1); 
        % This produces a column vector in which each of the 3 components 
        % is chosen from a normal distribution with a mean of 0 and 
        % standard deviation of 1.
        vector = KuhnLength*(RandomNumbers/norm(RandomNumbers));
        % This vector is of length KuhnLength and actually is pointing in a
        % totally random direction in 3D (i.e. the above expression is the 
        % the same as choosing the end of the vector from a uniform
        % distribution over a sphere of radius KuhnLength).
        % Interesting, if you have a vector in 3D whose components 
        % have been randomly chosen from the normal distribution mentioned above,
        % you can divide the vector by its magnitude in order to obtain 
        % a vector of unit length pointing in a totally random direction.
                
        CurrentCoordinates = CurrentCoordinates + vector;
        % This updates the current coordinates.
        CoordinateMatrix(:,u+2) = CurrentCoordinates;
        % This keeps track of the coordinates of each segment.        
        CDFC = norm(CurrentCoordinates);           
        u = u + 1; 
        if u == n
            GeneCoordinates = CurrentCoordinates;
            % This records the coordinates of the gene of interest.
        end
                
    end
 
    if (CDFC < R) && (CDFC > R - PeripheryThickness)...
         && (CurrentCoordinates(3,1) > -R + h) 
        GeneToSPBDistance = norm(GeneCoordinates - SPBCoordinates);
        Bin = ceil(GeneToSPBDistance*NumberOfBins/MaximumDistance);
        ProbabilityDistribution(1,Bin) = ProbabilityDistribution(1,Bin) + 1;
        SuccessfulConformation = CoordinateMatrix;
        % The "if statement" determines whether the three necessary
        % conditions have been met for the last segment. If they have,
        % then the distance between the gene and the SPB is calculated; 
        % then we determine which bin of the probability distribution this 
        % distance corresponds to, and we add a tally of 1 to that bin.               
    end
    
end

TimeInHours = toc/3600
NumberOfSuccessfulConformations = sum(ProbabilityDistribution)
% These tell the user how long the simulation took and how many chromosome
% conformations satisfied the three necessary conditions.

PDF = ProbabilityDistribution/(sum(ProbabilityDistribution)*dx); 
XDistances = (dx/2):dx:(MaximumDistance-dx/2);
figure(1)
hold on
plot(XDistances,PDF,'-dr','LineWidth',3)
xlabel('Gene to SPB Distance (micrometers)')
ylabel('PDF (1/micrometer)')
title('Probability Density Function for the Distance from the Gene to the SPB')
% The lines above calculate the probability density functions (PDF) by
% dividing the probability distribution by the width of each bin. The code
% then plots the PDF.

figure(2)
% The following code draws (in 3D!) the last successful conformation.
hold on
plot3(SuccessfulConformation(1,:),SuccessfulConformation(2,:),...
    SuccessfulConformation(3,:),'-m','LineWidth',3)
xlabel('x-coordinates (micrometers)')
ylabel('y-coordinates (micrometers)')
zlabel('z-coordinates (micrometers)')
xlim([-R R])
ylim([-R R])
zlim([-R R])
plot3(0,0,R,'o','MarkerEdgeColor','y','MarkerFaceColor','r','MarkerSize',20) 
plot3(SuccessfulConformation(1,n+1),SuccessfulConformation(2,n+1),...
    SuccessfulConformation(3,n+1),...
    'o','MarkerEdgeColor','y','MarkerFaceColor','g','MarkerSize',20) 
plot3(SuccessfulConformation(1,N+1),SuccessfulConformation(2,N+1),...
    SuccessfulConformation(3,N+1),...
    'o','MarkerEdgeColor','y','MarkerFaceColor','k','MarkerSize',20) 
[x,y,z] = sphere(50);
x = x*R;
y = y*R;
z = z*R;
lightGrey = 0.8*[1 1 1]; 
surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey)
X = [0 0 R -R];
Y = [R -R 0 0];
Z = (-R+h)*ones(4,4);
surf(X,Y,Z)
view(60,15)
axis equal
legend('Chromosome Arm','Spindle Pole Body','Gene of Interest',...
    'Telomere','Nuclear Envelope','Plane Indicating the Top of the Nucleolus')
FigHandle = figure(2);
  set(FigHandle, 'Position', [100, 100, 1100, 600]);
% Get different views of the conformation by using the rotate tool in the
% figure window. Also, compare the simulated conformation to an
% experimentally-determined conformation shown in Figure 5 of
% http://www.nature.com/nature/journal/v465/n7296/pdf/nature08973.pdf




