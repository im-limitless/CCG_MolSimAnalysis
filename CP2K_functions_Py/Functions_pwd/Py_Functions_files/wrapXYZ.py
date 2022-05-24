import numpy as np
    
def wrapXYZ(XYZ_input = None,Vec = None): 
    XYZ_wrap = XYZ_input
    for i in np.arange(1,XYZ_wrap.shape[1-1]+1).reshape(-1):
        for j in np.arange(1,XYZ_wrap.shape[2-1]+1).reshape(-1):
            if XYZ_wrap(i,j,1) < 0:
                XYZ_wrap[i,j,1] = XYZ_wrap(i,j,1) + Vec(1)
            else:
                if XYZ_wrap(i,j,1) > Vec(1):
                    XYZ_wrap[i,j,1] = XYZ_wrap(i,j,1) - Vec(1)
            if XYZ_wrap(i,j,2) < 0:
                XYZ_wrap[i,j,2] = XYZ_wrap(i,j,2) + Vec(2)
            else:
                if XYZ_wrap(i,j,2) > Vec(2):
                    XYZ_wrap[i,j,2] = XYZ_wrap(i,j,2) - Vec(2)
    
    return XYZ_wrap
    # # # attempt using logical indices has bug wrt 3d indexing, used loops
# above instead.
# # # # if any(any((XYZ_input(:,:,1) < 0),2))
# # #     XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) < 0) = XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) < 0) + Vec(1);
# # # #     disp('hello');
# # # # end
# # #
# # # # if any(any((XYZ_input(:,:,1) > Vec(1)),2))
# # # #     XYZ_input(XYZ_input(1,:,1) > Vec(1)) =  XYZ_input(XYZ_input(1,:,1)*ones(size(XYZ_input)) > Vec(1)) -  Vec(1);
# # # # XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) > Vec(1)) = XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) > Vec(1)) - Vec(1);
# # # #     disp('hello');
# # # # end
# # #
# # # # if any(any((XYZ_input(:,:,2) < 0),2))
# # # #     XYZ_input(XYZ_input(1,:,2) < 0) = XYZ_input(XYZ_input(1,:,2)*ones(size(XYZ_input)) < 0) + Vec(2);
# # # XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) < 0) = XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) < 0) + Vec(2);
# # # #     disp('hello');
# # # # end
# # #
# # # # if any(any((XYZ_input(:,:,2) > Vec(2)),2))
# # # #     XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) > Vec(2)) = XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) > Vec(2)) - Vec(2);
# # # #     disp('hello');
# # # # end
# # #
# # # XYZ_wrap = XYZ_input;
    
    # return
    return XYZ_wrap