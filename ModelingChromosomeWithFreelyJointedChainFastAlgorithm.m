% Modeling a Chromosome with a Freely Jointed Chain
% Using the FAST ALGORITHM for Simulating A Confined and/or Tethered Chromosomes
% Gabe Bronk, 2016
% gbronk@brandeis.edu
% Kondev Group, Brandeis University 


% This MATLAB code simulates a random conformation of the arm of a chromosome 
% in a yeast cell's nucleus. This script is intended for any researcher who 
% wishes to model a chromosome in a simple manner and wants to reduce the
% computation time compared to doing a regular random walk chromosome simulation.
% This code uses the FAST ALGORITHM: This algorithm speeds up 
% the simulation of a freely-jointed chain chromosome that is confined to a 
% certain region and/or attached to certain locations. This algorithm was 
% developed by Bronk for the paper "Chromosome-refolding model of mating-type
% switching in yeast" published in PNAS in 2016 by Avaroglu, Bronk, et al.
% The algorithm is described in the paper's Supplementary Information under
% the heading "Polymer Model Simulation Details."

% Prior to using this script, first look over the MATLAB script entitled
% "IntroductionToModelingChromosomesWithFreelyJointedChain.m" which is
% simpler and will introduce you to the basics. After looking at that
% script, come back to this script and look at the comments that have a 
%%% triple percentage symbol "%%%" next to them. This indicates the comments 
%%% specifically pertaining to the FAST ALGORITHM.

% This MATLAB code simulates a random conformation of the arm of a chromosome 
% in a yeast cell's nucleus. The chromosome is modelled as many straight
% segments, with each segment oriented in a random direction. It's like 
% drawing a zig-zagging line (where each segment is the same length), but
% the line is zig-zagging in 3D. The nucleus is modelled as a sphere, 
% and the nucleolus is modelled as a spherical cap.

%%% Summary of the FAST ALGORITHM: 
%%% When performing random walks to generate chromosome conformations it is 
%%% exceedingly rare to generate a conformation that fulfills
%%% all of the constraints described in "Chromosome-refolding model of mating-type
%%% switching in yeast" in the section entitled Materials and Methods, 
%%% Computational and Theoretical Methods, Freely Jointed Chain Simulation.
%%% Examples of these constraints include ensuring that the chromosome does not
%%% leave the nucleus or enter the nucleolus region, or ensuring that a
%%% particular chromosomal locus is attached to a particular point in
%%% space or attached to any location on the nuclear membrane.  
%%% To decrease the otherwise prohibitive simulation time, we developed an
%%% algorithm in which we alter conformations to fit the constraints.
%%% An arm of a chromosome is simulated as a random walk, starting from the
%%% centromere (located 50 nm below the north pole of the nuclear
%%% sphere). During each step of the random walk, the Kuhn segment is
%%% reoriented if it exits the nucleus or enters the nucleolus and is also
%%% reoriented so that particular chromosomal loci are positioned according to 
%%% your desired constraints (There are places in this code where you can input
%%% your desired constraints). This reorientation, however, would result 
%%% in the generation of some allowed conformations more frequently than 
%%% others, in violation of thermodynamics which requires that all freely  
%%% jointed chain conformations that satisfy the constraints be equally likely. 
%%% We therefore calculate the resulting
%%% increase in the probability of each conformation and assign to each
%%% conformation a weight to counteract the increase in the probability.

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
NumberOfTrials = 10^4; 
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

%%% PeripheryThickness = 0.05; % units in micrometers
%%% This note is in regards to any membrane-bound locus (the telomere,
%%% for instance, is bound to the membrane in yeast).
%%%% "PeripheryThickness" has been removed from this script because the
%%%% this fast algorithm actually places the membrane-bound locus exactly
%%%% on the membrane (i.e. at a distance R from the center).
MicrotubuleLength = 0.05; % units in micrometers
% This is the length of the microtubule that connects the chromosome's
% centromere to the SPB.
NumberOfBins = 50;
% Number of bins of the probability distribution.

%%% The following matrices take inputs of what parts of the chromosome
%%% should be attached to certain locations. Here's an example:
%%% ChromatinToAttach = [0.020 0.048]; 
%%% PlacesToAttachTo = [[0; 0; R] [0; 0; 0]];
%%% In ChromatinToAttach, list (in units of megabase pairs from the
%%% centromere) which parts of the chromosome you would like to be attached
%%% to places. You must write these in ascending order - i.e. don't write
%%% ChromatinToAttach = [0.048 0.020]; 
%%% In PlacesToAttachTo, list the coordinates of the places that you want those
%%% parts of the chromosome to attach to. In the example above, the part of
%%% the chromosome 0.020 megabase pairs from the centromere will be attached
%%% to the coordinate (0,0,R). And the the part of
%%% the chromosome 0.048 megabase pairs from the centromere will be attached
%%% to the coordinate (0,0,0).
%%% If you don't want to attach any loci to anything, then write:
%%% ChromatinToAttach = []; 
%%% PlacesToAttachTo = [];
ChromatinToAttach = [0.070]; %This number is just an example - it has no biological backing.
PlacesToAttachTo = [[R; 0; 0]];

%%% The following matrix takes inputs of what parts of the chromosome
%%% should be attached to anywhere on the nuclear membrane. Here's an example:
%%% MembraneBoundChromatin = [0.070 0.100 0.112 0.124]; 
%%% This lists (in units of megabase pairs from the
%%% centromere) which parts of the chromosome you would like to be attached
%%% to anywhere on the nuclear membrane.
%%% If you don't want any membrane-bound loci, then write:
%%% MembraneBoundChromatin = []; 
MembraneBoundChromatin = [0.124]; % This example makes the telomere membrane-bound
                                  % since we also set LengthFromCentromereToTelomere to equal 0.124.

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
% Number of Kuhn segments for the whole chromosome arm.
n = round((ChromosomeCompaction*LengthFromCentromereToGene)/((KuhnLength)^2));
% Number of Kuhn segments between the centromere and the gene of interest.

CommandArray = ones(1,N);
for i = 1:length(ChromatinToAttach)
AttachedSegmentNumber(i) = round((ChromosomeCompaction*ChromatinToAttach(i))/((KuhnLength)^2));
CommandArray(AttachedSegmentNumber(i) - 1) = 2;
CommandArray(AttachedSegmentNumber(i)) = 3;
end
for i = 1:length(MembraneBoundChromatin)
MembraneBoundSegmentNumber(i) = round((ChromosomeCompaction*MembraneBoundChromatin(i))/((KuhnLength)^2));
CommandArray(MembraneBoundSegmentNumber(i)) = 4;
end


MaximumDistance = 2*R; 
% The maximum distance to include in the probability distribution
dx = MaximumDistance/NumberOfBins;
% The width of each bin of the probability distribution.
ProbabilityDistribution = zeros(1,NumberOfBins);
% This starts the probabiity distribution at all zeros, and we will
% increase the probabilities in the bins as the simulation progresses.
InterGeneProbabilityDistribution = zeros(1,NumberOfBins);


hroots = roots([(pi/3) (-pi*R) 0 (Nucleolus*4/3*pi*R^3)]);
hindex = find(hroots > 0 & hroots < (2*R));
h = hroots(hindex);
% h is the height of the nucleolus (in micrometers). The lines above do the
% following calculation: We represent the nucleolus as a spherical cap,
% and we want to set this spherical cap to a particular volume:
% e.g. 1/5 of the total nuclear volume. What must the height of the spherical
% cap be in order for the spherical cap to have that particular volume? The
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
%%% We can also, for whatever reason, impose other constraints, namely
%%% attaching loci to certain locations.

% The number of times that
% we try to simulate the chromosome arm is "NumberOfTrials". However, most
% of these attempts will fail. This is because most random conformations do not 
% have their telomeres close to the nuclear envelope and/or they protrude 
% out of the nuclear sphere and/or penetrate into the nucleolus.  
% Note that in reality, some protrusion into
% the nucleolus is possible, but our simplified model represents the
% nucleolus as impenetrable).

NumberOfSuccessfulConformations = 0;
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
        
    SegmentWeights = zeros(1, N);
    %%% This keeps track of the weight of each kuhn segment. The weights
    %%% are for the weighting scheme that allows the simulation to run
    %%% faster.
    AttachmentTally = 0;
    StopTheConformation = 0;
    %%% This line simply allows the while loop below to start.
    
    while (u < N) && StopTheConformation == 0
        %%% The "while loop" tells MATLAB to keep adding segments until
        %%% it has added N segments. Also, it tells MATLAB to stop adding 
        %%% segments if the previous semgent violated one of the necessary 
        %%% conditions (e.g. the previous segment has protruded out of the nuclear
        %%% sphere or into the nucleolus, or was too far from a point of attachment).
        %%% If a segment does such violation,
        %%% the while loop does not proceed, the data is not recorded 
        %%% (i.e. we throw away all the segments of this conformation), and
        %%% we go to the next trial of the "for loop" in order to start a
        %%% brand new conformation.
        %%%
        %%% StopTheConformation = 0 means that no violations of the necessary conditions
        %%% have occured. StopTheConformation == 1 means that violations have
        %%% occured.
        
        %%% All lines contained in the following "if statement" assign
        %%% weights to the kuhn segments. A particular weight is assigned 
        %%% to a kuhn segment depending on the kuhn segment's location and
        %%% depending on what we are forcing the kuhn segment to do (i.e.
        %%% are we attaching the end of a kuhn segment to a particular spot in space, 
        %%% or are we forcing the end of the kuhn segment to reside
        %%% anywhere on the nuclear membrane? Are we stopping the kuhn
        %%% segment from leaving the nucleus or from entering the
        %%% nucleolus?). The weights were derived for particular
        %%% situations using algebra and geometry. Make sure not to change
        %%% the weights since that would make the simulation incorrect.
        %%% 
        %%% The lines below also determine the orientation of the next kuhn
        %%% segment.
             
        
            if CommandArray(u+1) == 1 
                if (CDFC > R - KuhnLength) 
                    SegmentWeights(u+1) =...
                        (KuhnLength + (R^2 - KuhnLength^2 - CDFC^2)/(2*CDFC))/(2*KuhnLength);
                    CDFC = R + 1; % This allows the following while loop to start.
                    while CDFC > R  %%% This loop is for reorienting Kuhn segments that exit the nucleus.
                        RandomNumbers = randn(3,1);
                        vector = KuhnLength/norm(RandomNumbers)*RandomNumbers;
                        PreliminaryCoordinates = CurrentCoordinates + vector;
                        CDFC = norm(PreliminaryCoordinates);
                    end
                    if (CurrentCoordinates(3,1) < -R + h)
                        StopTheConformation = 1;
                    end
                else if (CDFC <= R - KuhnLength) && (CurrentCoordinates(3,1) < -R + h + KuhnLength)
                        SegmentWeights(u+1) =...
                            (KuhnLength + R + CurrentCoordinates(3,1) - h)/(2*KuhnLength);
                        z = -R; % This allows the following while loop to start.
                        while (z < -R + h) %%% This loop is for reorienting Kuhn segments that protrude into the nucleolus.
                            RandomNumbers = randn(3,1);
                            vector = KuhnLength/norm(RandomNumbers)*RandomNumbers;
                            PreliminaryCoordinates = CurrentCoordinates + vector;
                            z = PreliminaryCoordinates(3,1);
                        end
                else
                        SegmentWeights(u+1) = 1;
                        RandomNumbers = randn(3,1);
                        vector = KuhnLength/norm(RandomNumbers)*RandomNumbers;
                        PreliminaryCoordinates = CurrentCoordinates + vector;
                end; end
      
            %%%------------------------------------------------------------    
            %%% The 32 lines after this comment are for forcing 
            %%% a particular chromosomal locus of interest 
            %%% to be attached to a particular point in space. 
            %%% However, not every chromosome conformation will be able to be constrained 
            %%% to attach to the desired point in space. What we do is we perform the random
            %%% walk that constructs the chromosome, and when we get to the 
            %%% Kuhn segment that comes 2 Kuhn segments before the locus of interest,
            %%% we ask: Is this kuhn segment within 2 Kuhn lengths of the 
            %%% desired point in space? 
                    %%% If the answer is "No": 
                    %%% We throw away the entire chromosome conformation 
                    %%% because there is no way to get the locus of interest to 
                    %%% be located at the desired point in space. 
                    %%% If the answer is "Yes": 
                    %%% We are then able to make sure to orient
                    %%% the next two Kuhn segments so that the locus of interest
                    %%% lands precisely on the desired point in space. Note that
                    %%% the Kuhn segment before the locus of interest will still be
                    %%% oriented randomly, but we constrain it such that its end
                    %%% is at any location that is exactly one Kuhn length from 
                    %%% the desired point in space.
            else if CommandArray(u+1) == 2 
                AttachmentTally = AttachmentTally + 1;
                d = norm(CurrentCoordinates - PlacesToAttachTo(:,AttachmentTally));
                if d >= 2*KuhnLength
                    StopTheConformation = 1;
                else if d < 2*KuhnLength 
                    a = sqrt(KuhnLength^2 - (d^2)/4);
                    F = PlacesToAttachTo(1,AttachmentTally) - CurrentCoordinates(1,1);
                    G = PlacesToAttachTo(2,AttachmentTally) - CurrentCoordinates(2,1);
                    H = PlacesToAttachTo(3,AttachmentTally) - CurrentCoordinates(3,1);
                    B = ( ((F^2/G + G)/(-H))^2 + 1 + F^2/G^2 )^(-0.5);
                    BasisVector1 = [sqrt(1 - F^2/(G^2 + F^2)); -F/sqrt(G^2 + F^2); 0];
                    BasisVector2 = [F/G*B; B; (F^2/G + G)/(-H)*B];
                    if (BasisVector1'*BasisVector2 > 10^-10) || (BasisVector1'*BasisVector2 < -1*10^-10)
                        BasisVector1 = [sqrt(1 - F^2/(G^2 + F^2)); F/sqrt(G^2 + F^2); 0];
                    end
                    RandomAngle = 2*pi*rand(1,1);
                    CurrentCoordinates = CurrentCoordinates + [F/2; G/2; H/2] +...
                        a*cos(RandomAngle)*BasisVector1 + a*sin(RandomAngle)*BasisVector2;
                    CDFC = norm(CurrentCoordinates);
                        
                    if (CurrentCoordinates(3,1) >= -R + h) && (CDFC < R)
                        SegmentWeights(u+1) = sqrt(KuhnLength^2 - d^2/4)/KuhnLength;
                    else
                        StopTheConformation = 1;
                    end
                end; end
                    
        else if CommandArray(u+1) == 3
            SegmentWeights(u+1) = 1;
            PreliminaryCoordinates = PlacesToAttachTo(:,AttachmentTally);
        %%%----------------------------------------------------------------
        else if CommandArray(u+1) == 4        
            if (CDFC <= R - KuhnLength)
                StopTheConformation = 1;
            else %%% This "else" statment is for attaching a particular chromosomal locus
                 %%% anywhere on the nuclear membrane, provided that 
                 %%% the previous Kuhn segment was already close enough to 
                 %%% the nuclear membrane that the current Kuhn segment can 
                 %%% be reoriented to attach to the nuclear membrane.
                d = R - CDFC;
                a = sqrt(KuhnLength^2 - ((2*R*d - d^2 - KuhnLength^2)/(2*R - 2*d))^2);
                F = CurrentCoordinates(1,1);
                G = CurrentCoordinates(2,1);
                H = CurrentCoordinates(3,1);
                B = ( ((F^2/G + G)/(-H))^2 + 1 + F^2/G^2 )^(-0.5);
                BasisVector1 = [sqrt(1 - F^2/(G^2 + F^2)); -F/sqrt(G^2 + F^2); 0];
                BasisVector2 = [F/G*B; B; (F^2/G + G)/(-H)*B];
                if (BasisVector1'*BasisVector2 > 10^-10) || (BasisVector1'*BasisVector2 < -1*10^-10)
                    BasisVector1 = [sqrt(1 - F^2/(G^2 + F^2)); F/sqrt(G^2 + F^2); 0];
                end
                RandomAngle = 2*pi*rand(1,1);
                PreliminaryCoordinates = CurrentCoordinates + CurrentCoordinates/norm(CurrentCoordinates)*(2*R*d-d^2-KuhnLength^2)/(2*R - 2*d) +...
                    a*cos(RandomAngle)*BasisVector1 + a*sin(RandomAngle)*BasisVector2;
                if (CurrentCoordinates(3,1) >= -R + h)
                    SegmentWeights(u+1) = a/KuhnLength;
                else
                    StopTheConformation = 1;
                end
            end
                
        end; end; end; end
     
        %%% Now that the kuhn segment's weight has been assigned (or the conformation has
        %%% been deemed in violation of the necessary conditions), we
        %%% record the coordinates of the Kuhn segment below.
        if StopTheConformation == 0
            CurrentCoordinates = PreliminaryCoordinates;
            CDFC = norm(CurrentCoordinates);
            CoordinateMatrix(:,u+2) = CurrentCoordinates;
            if u == n
                GeneCoordinates = CurrentCoordinates;
                % This records the coordinates of the gene of interest.
            end
            u = u + 1;
        end
        
    end         
              
 
    if StopTheConformation == 0
        %%% This "if statement" determines whether no kuhn segments have violated
        %%% the necessary conditions. If no violations have occured,
        %%% then the distance between the gene and the SPB is calculated; 
        %%% then we determine which bin of the probability distribution this 
        %%% distance corresponds to, and we add the weight of that conformation
        %%% to that bin. 
        GeneToSPBDistance = norm(GeneCoordinates - SPBCoordinates);
        Bin = ceil(GeneToSPBDistance*NumberOfBins/MaximumDistance);
        TotalConformationWeight = prod(SegmentWeights);
        ProbabilityDistribution(1,Bin) =...
            ProbabilityDistribution(1,Bin) + TotalConformationWeight;
        SuccessfulConformation = CoordinateMatrix;
        SegmentWeightsOfSuccessfulConformation = SegmentWeights;
        %%% The above line can be helpful for debugging.
        NumberOfSuccessfulConformations = NumberOfSuccessfulConformations + 1;
        %%% The above line keeps track of the number of successful
        %%% conformations simulated. However, due to the weighting scheme of
        %%% the fast algorithm, this number must be taken with a grain of
        %%% salt: 1000 successful conformations with the fast algorithm does
        %%% not give as smooth a distribution as 1000 successful
        %%% conformations with the slow algorithm.
        
        %%% The following "if statement" records the coordinates of a gene
        %%% from one trial: call the coordinates (x, y, z). 
        %%% Then the gene's coordinates in the next successful trial
        %%% (call these coordinate (a,b,c)) are found, and the code calculates the
        %%% distance between the genes in the two trials (i.e. it finds the
        %%% distance between (x,y,z) and (a,b,c).
        %%%  This could be useful in diploid
        %%% drosophila nuclei in which two GFP spots are seen, each GFP spot 
        %%% marking each of two copies of the same gene.
        if mod(NumberOfSuccessfulConformations,2) == 1
            SavedGeneCoordinates = GeneCoordinates;
            SavedTotalConformationWeight = TotalConformationWeight;
        else
            GeneToGeneDistance = norm(GeneCoordinates - SavedGeneCoordinates);
            Bin = ceil(GeneToSPBDistance*NumberOfBins/MaximumDistance);
            WeightForBothConformations =...
                TotalConformationWeight*SavedTotalConformationWeight;
            InterGeneProbabilityDistribution(1,Bin) =...
                InterGeneProbabilityDistribution(1,Bin) + WeightForBothConformations;
        end
        
    end
    
end

TimeInHours = toc/3600
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

InterGenePDF = InterGeneProbabilityDistribution/(sum(InterGeneProbabilityDistribution)*dx); 
figure(3)
hold on
plot(XDistances,InterGenePDF,'-dg','LineWidth',3)
xlabel('Gene to Gene Distance (micrometers)')
ylabel('PDF (1/micrometer)')
title('Probability Density Function for the Distance Between the Two Copies of the Gene')
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




