function [XYZ_wrap] = wrapXYZ(XYZ_input, Vec)

XYZ_wrap = XYZ_input;

for i = 1:size(XYZ_wrap,1)
    for j = 1:size(XYZ_wrap,2)
        if XYZ_wrap(i,j,1) < 0-0.1*Vec(1)
            XYZ_wrap(i,j,1) = XYZ_wrap(i,j,1) + Vec(1);
        elseif XYZ_wrap(i,j,1) > Vec(1)+0.05*Vec(1)
            XYZ_wrap(i,j,1) = XYZ_wrap(i,j,1) - Vec(1);
        end
        
        if XYZ_wrap(i,j,2) < 0-0.1*Vec(2)
            XYZ_wrap(i,j,2) = XYZ_wrap(i,j,2) + Vec(2);
        elseif XYZ_wrap(i,j,2) > Vec(2)+0.05*Vec(2)
            XYZ_wrap(i,j,2) = XYZ_wrap(i,j,2) - Vec(2);
        end
    end
end

return
% % % attempt using logical indices has bug wrt 3d indexing, used loops
% above instead.
% % % % if any(any((XYZ_input(:,:,1) < 0),2))
% % %     XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) < 0) = XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) < 0) + Vec(1);
% % % %     disp('hello');
% % % % end
% % % 
% % % % if any(any((XYZ_input(:,:,1) > Vec(1)),2))
% % % %     XYZ_input(XYZ_input(1,:,1) > Vec(1)) =  XYZ_input(XYZ_input(1,:,1)*ones(size(XYZ_input)) > Vec(1)) -  Vec(1);
% % % % XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) > Vec(1)) = XYZ_input(XYZ_input(:,:,1).*ones(size(XYZ_input)) > Vec(1)) - Vec(1);
% % % %     disp('hello');
% % % % end
% % % 
% % % % if any(any((XYZ_input(:,:,2) < 0),2))
% % % %     XYZ_input(XYZ_input(1,:,2) < 0) = XYZ_input(XYZ_input(1,:,2)*ones(size(XYZ_input)) < 0) + Vec(2);
% % % XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) < 0) = XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) < 0) + Vec(2);
% % % %     disp('hello');
% % % % end
% % % 
% % % % if any(any((XYZ_input(:,:,2) > Vec(2)),2))
% % % %     XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) > Vec(2)) = XYZ_input(XYZ_input(:,:,2).*ones(size(XYZ_input)) > Vec(2)) - Vec(2);
% % % %     disp('hello');
% % % % end
% % % 
% % % XYZ_wrap = XYZ_input;

% return
