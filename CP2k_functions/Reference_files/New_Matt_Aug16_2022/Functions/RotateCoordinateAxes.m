function [XYZR] = RotateCoordinateAxes(XYZ, ABC, theta, Axes, TF)
% theta is angle of rotation. Axes is the cartesian axis to rotate about x = 1, y =
% 2, z = 3. TF is true/fale, 1 for rotate 0 for do nothing.

% XYZR = zeros(size(XYZ));

if TF
    switch Axes
        case 1
            Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            XYZR = reshape(XYZ, [size(XYZ ,1)*size(XYZ,2) size(XYZ,3)])*Rx;
            XYZR = reshape(XYZR, [size(XYZ,1) size(XYZ,2) size(XYZ,3)]);
%             XYZR(:,:,1) = XYZ(:,:,1).*Rx(1,1)+XYZ(:,:,2).*Rx(1,2)+XYZ(:,:,3).*Rx(1,3);
%             XYZR(:,:,2) = XYZ(:,:,1).*Rx(2,1)+XYZ(:,:,2).*Rx(2,2)+XYZ(:,:,3).*Rx(2,3);
%             XYZR(:,:,3) = XYZ(:,:,1).*Rx(3,1)+XYZ(:,:,2).*Rx(3,2)+XYZ(:,:,3).*Rx(3,3);
%         case 2
%             Ry = [];
%             XYZR(:,:,1) = XYZ(:,:,1).*Ry(1,1)+XYZ(:,:,2).*Ry(1,2)+XYZ(:,:,3).*Ry(1,3);
%             XYZR(:,:,2) = XYZ(:,:,1).*Ry(2,1)+XYZ(:,:,2).*Ry(2,2)+XYZ(:,:,3).*Ry(2,3);
%             XYZR(:,:,3) = XYZ(:,:,1).*Ry(3,1)+XYZ(:,:,2).*Ry(3,2)+XYZ(:,:,3).*Ry(3,3);
%         case 3
%             Rz = [];
%             XYZR(:,:,1) = XYZ(:,:,1).*Rz(1,1)+XYZ(:,:,2).*Rz(1,2)+XYZ(:,:,3).*Rz(1,3);
%             XYZR(:,:,2) = XYZ(:,:,1).*Rz(2,1)+XYZ(:,:,2).*Rz(2,2)+XYZ(:,:,3).*Rz(2,3);
%             XYZR(:,:,3) = XYZ(:,:,1).*Rz(3,1)+XYZ(:,:,2).*Rz(3,2)+XYZ(:,:,3).*Rz(3,3);
%             
    end
else
    XYZR = XYZ;
    return
end