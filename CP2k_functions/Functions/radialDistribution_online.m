function gR = radialDistribution(switchVal,gR,coords,L)
    
        % Set three operation cases for this function
        initialize=0;
        sample=1;
        results=2;
    
        % Choose the operation method according to the providerd switchVal
        switch switchVal
            
            case initialize
                % Initialize a histogram to hold the radial distribution
                
                gR.count = 0;
                gR.range = [0 0.5*L]; % L is length of box
                gR.increment = L/200.0; % bin width
                gR.outFreq = 1000;
                gR.saveFileName = 'radialDist.dat';
                
            case sample
                % Loop over pairs and determine the distribution of distances
                nPart = size(coords,2);
                
                for partA = 1:(nPart-1)
                    for partB = (partA+1):nPart
                        % Calculate particle-particle distance
                        % Account for PBC (assuming 3D)                               
                        dr = coords(:,partA) - coords(:,partB); 
                        dr = distPBC3D(dr,L);
                        % Get the size of this distance vector
                        r = sqrt(sum(dot(dr,dr)));
                        
                        % Add to g(r) if r is in the right range [0 L/2]
                        if (r < 0.5*L)
                            gR = histogram(gR,r);
                        end
                    end
                end
                
            case results
                % The radial distribution function should be normalized.
                % First, just like any other hisotgram:
                gR.histo = gR.histo/(gR.count*gR.increment);
                
                % Now, each bin should be normalized according to its volume
                % since larger inter-particle distances are expected to contain 
                % more counts even in the ideal case. 
                % For more information see Frenkel & Smit chapter 4.4
                
                nBins = size(gR.values,2);
                nPart = size(coords,2); % (length(IndxO)+length(IndxH))
                rho = nPart/(L^3); % Density of the particles npart/(ABC(1)*ABC(2)*Pt)
                
                for bin = 1:nBins
                    rVal = gR.values(bin);
                    next_rVal = rVal+gR.increment;
                    % Calculate the volume of the bin
                    volBin = (4/3.0)*pi*(next_rVal^3 - rVal^3);
                    % Calculate the number of particles expected in this bin in
                    % the ideal case
                    nIdeal = volBin*rho;
                    % Normalize the bin
                    gR.histo(bin) = gR.histo(bin) / nIdeal;
                end
                
            otherwise
                % Wrong switch
                disp('radialDistribution : You have entered an illegal switch value');
        end
    
    end