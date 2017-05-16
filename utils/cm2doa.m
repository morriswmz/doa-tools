function doa = cm2doa(z, k, unit)
%CM2DOA Converts complex exponentials to doas.
switch lower(unit)
    case 'radian'
        doa = asin(angle(z) / k);
    case 'degree'
        doa = rad2deg(asin(angle(z) / k));
    case 'sin'
        doa = angle(z) / k;
    otherwise
        error('Unkown unit ''%s''.', unit);
end
end

