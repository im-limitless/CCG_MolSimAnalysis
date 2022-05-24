import numpy as np
    
def gaussfilt(t = None,z = None,sigma = None): 
    #Apply a Gaussian filter to a time series
#   Inputs: t = independent variable, z = data at points t, and
#       sigma = standard deviation of Gaussian filter to be applied.
#   Outputs: zfilt = filtered data.
    
    #   written by James Conder. Aug 22, 2013
#   Sep 04, 2014: Convolution for uniformly spaced time time vector (faster)
#   Mar 20, 2018: Damped edge effect of conv (hat tip to Aaron Close)
    
    n = len(z)
    
    a = 1 / (np.sqrt(2 * pi) * sigma)
    
    sigma2 = sigma * sigma
    # check for uniform spacing
# if so, use convolution. if not use numerical integration
    uniform = False
    dt = np.diff(t)
    dt = dt(1)
    ddiff = np.amax(np.abs(np.diff(np.diff(t))))
    if ddiff / dt < 0.0001:
        uniform = True
    
    if uniform:
        filter = dt * a * np.exp(- 0.5 * ((t - mean(t)) ** 2) / (sigma2))
        i = filter < dt * a * 1e-06
        filter[i] = []
        zfilt = conv(z,filter,'same')
        onesToFilt = np.ones((z.shape,z.shape))
        onesFilt = conv(onesToFilt,filter,'same')
        zfilt = zfilt / onesFilt
    else:
        ### get distances between points for proper weighting
        w = 0 * t
        w[np.arange[2,end() - 1+1]] = 0.5 * (t(np.arange(3,end()+1)) - t(np.arange(1,end() - 2+1)))
        w[1] = t(2) - t(1)
        w[end()] = t(end()) - t(end() - 1)
        ### check if sigma smaller than data spacing
        iw = find(w > 2 * sigma,1)
        if not len(iw)==0 :
            print('WARNING: sigma smaller than half node spacing')
            print('May lead to unstable result')
            iw = w > 2.5 * sigma
            w[iw] = 2.5 * sigma
            # this correction leaves some residual for spacing between 2-3sigma.
# otherwise ok.
# In general, using a Gaussian filter with sigma less than spacing is
# a bad idea anyway...
        ### loop over points
        zfilt = 0 * z
        for i in np.arange(1,n+1).reshape(-1):
            filter = a * np.exp(- 0.5 * ((t - t(i)) ** 2) / (sigma2))
            zfilt[i] = sum(np.multiply(np.multiply(w,z),filter))
        ### clean-up edges - mirror data for correction
        ss = 2.4 * sigma
        # left edge
        tedge = np.amin(t)
        iedge = find(t < tedge + ss)
        nedge = len(iedge)
        for i in np.arange(1,nedge+1).reshape(-1):
            dist = t(iedge(i)) - tedge
            include = find(t > t(iedge(i)) + dist)
            filter = a * np.exp(- 0.5 * ((t(include) - t(iedge(i))) ** 2) / (sigma2))
            zfilt[iedge[i]] = zfilt(iedge(i)) + sum(np.multiply(np.multiply(w(include),filter),z(include)))
        # right edge
        tedge = np.amax(t)
        iedge = find(t > tedge - ss)
        nedge = len(iedge)
        for i in np.arange(1,nedge+1).reshape(-1):
            dist = tedge - t(iedge(i))
            include = find(t < t(iedge(i)) - dist)
            filter = a * np.exp(- 0.5 * ((t(include) - t(iedge(i))) ** 2) / (sigma2))
            zfilt[iedge[i]] = zfilt(iedge(i)) + sum(np.multiply(np.multiply(w(include),filter),z(include)))
    
    return zfilt
    
    return zfilt