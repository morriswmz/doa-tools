function s = seconds2str(t)
%SECONDS2STR Convert seconds to nicer looking string.
t = double(t);
if t < 1e-3
    s = sprintf('%3.2f ns', t*1e6);
elseif t < 1
    s = sprintf('%3.2f ms', t*1e3);
elseif t < 60
    s = sprintf('%2.2f s', t);
else
    if t > 31536000
        year = floor(t/31536000);
        t = mod(t, 31536000);
    else
        year = 0;
    end
    if t > 86400
        day = floor(t/86400);
        t = mod(t, 86400);
    else
        day = 0;
    end
    if t > 3600
        hour = floor(t/3600);
        t = mod(t, 3600);
    else
        hour = 0;
    end
    minute = floor(t/60);
    second = floor(mod(t, 60));
    s = sprintf('%d m %d s', minute, second);
    if hour > 0
        s = [sprintf('%d h ', hour) s];
    end
    if day > 0
        s = [sprintf('%d d ', day) s];
    end
    if year > 0
        s = [sprintf('%d y ', day) s];
    end
end

