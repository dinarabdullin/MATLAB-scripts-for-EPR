function [dir] = randomDistrOnSphere(size)

%     size = 1000000;
    xi = acos(2 * rand(size,1) - 1);
    phi = 2 * pi * rand(size,1); 
%     figure;
%     subplot(1,2,1);
%     hist(xi *180/pi,180);
%     xlabel('Xi (deg)');
%     ylabel('Probability');
%     subplot(1,2,2);
%     hist(phi * 180/pi,360);
%     xlabel('Phi (deg)');
%     ylabel('Probability');
    dir = zeros(3,size);
    for i = 1:size 
        dir(:,i) = spherical2Cartesian(1.0, xi(i), phi(i));
    end

end