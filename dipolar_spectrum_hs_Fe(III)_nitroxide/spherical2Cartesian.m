function [v] = spherical2Cartesian(r, theta, phi)

    v = zeros(3,1);
    v(1,1) = r * cos(phi) * sin(theta);
    v(2,1) = r * sin(phi) * sin(theta);
    v(3,1) = r * cos(theta);

end

