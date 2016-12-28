function b = ae2broad(azimuth, elevation)
%AE2BROAD Converts (azimuth, elevation) angles to broadside angles.
b = asin(sin(azimuth).*cos(elevation));
end

