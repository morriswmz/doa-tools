function x = randccsn(varargin)
%RANDCCSN Generates complex circularly-symmetric normal random numbers.
%Syntax follows randn
x = sqrt(1/2) * (randn(varargin{:}) + 1j*randn(varargin{:}));
end

