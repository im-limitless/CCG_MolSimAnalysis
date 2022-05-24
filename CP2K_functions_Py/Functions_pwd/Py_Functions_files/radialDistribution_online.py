import numpy as np
    
def radialDistribution(switchVal = None,gR = None,coords = None,L = None): 
    # Set three operation cases for this function
    initialize = 0
    sample = 1
    results = 2
    # Choose the operation method according to the providerd switchVal
    if initialize == switchVal:
        # Initialize a histogram to hold the radial distribution
        gR.count = 0
        gR.range = np.array([0,0.5 * L])
        gR.increment = L / 200.0
        gR.outFreq = 1000
        gR.saveFileName = 'radialDist.dat'
    else:
        if sample == switchVal:
            # Loop over pairs and determine the distribution of distances
            nPart = coords.shape[2-1]
            for partA in np.arange(1,(nPart - 1)+1).reshape(-1):
                for partB in np.arange((partA + 1),nPart+1).reshape(-1):
                    # Calculate particle-particle distance
# Account for PBC (assuming 3D)
                    dr = coords(:,partA) - coords(:,partB)
                    dr = distPBC3D(dr,L)
                    # Get the size of this distance vector
                    r = np.sqrt(sum(np.dot(dr,dr)))
                    # Add to g(r) if r is in the right range [0 L/2]
                    if (r < 0.5 * L):
                        gR = histogram(gR,r)
        else:
            if results == switchVal:
                # The radial distribution function should be normalized.
# First, just like any other hisotgram:
                gR.histo = gR.histo / (gR.count * gR.increment)
                # Now, each bin should be normalized according to its volume
# since larger inter-particle distances are expected to contain
# more counts even in the ideal case.
# For more information see Frenkel & Smit chapter 4.4
                nBins = gR.values.shape[2-1]
                nPart = coords.shape[2-1]
                rho = nPart / (L ** 3)
                for bin in np.arange(1,nBins+1).reshape(-1):
                    rVal = gR.values(bin)
                    next_rVal = rVal + gR.increment
                    # Calculate the volume of the bin
                    volBin = (4 / 3.0) * pi * (next_rVal ** 3 - rVal ** 3)
                    # Calculate the number of particles expected in this bin in
# the ideal case
                    nIdeal = volBin * rho
                    # Normalize the bin
                    gR.histo[bin] = gR.histo(bin) / nIdeal
            else:
                # Wrong switch
                print('radialDistribution : You have entered an illegal switch value')
    
    return gR
    
    return gR