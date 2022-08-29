function [v,h] = GETFOURIERCOEFF( h, x_t )
%GETFOURIERCOEFF.m is a function that takes the inputs of a function in the
%time domain and the harmonics of interest, and returns the fourier
%coefficients pertaining to the harmonics in the vector h of thefunction
%x_t.
% USAGE:
%   [v,h] = GETFOURIERCOEFF(h,x_t)
% INPUTS:
%   h		: A vector containing the harmonics pertaining to the fourier
%      		  coefficients that will be returned in vector v. 0
%      		  will be pushed to the top if it's not already there.
%   x_t		: A function in the time domain of which the fourier coefficients
%        	  will be found.
% OUTPUTS:
%   v		: A vector containing the fourier coefficients pertaining to the
%      		  harmonics found in vector h
%   h		: List of harmonics (unchanged if 0 was not
%   		  there/was in the first position)
% EXAMPLE:
% t = linspace(0, 2*pi, 129); t(end) = [];
% y = cos(t') + sin(2*t') + 4*sin(3*t');
% yharmonics = GETFOURIERCOEFF([0; 1; 2; 3; 4], y);
% % The above will return a 9x1 vector of harmonics as:
% % [y0; yc1; ys1; yc2; ys2; yc3; ys3; yc4; ys4]
% % with y0 being 0th harmonic; yc's being cosine harmonics and ys's being sine harmonics

%     % IMPLEMENTATION REQUIRING THE PRESENCE OF THE ZERO HARMONIC ALWAYS
%     n = length(h)-1;
%     nd = size(x_t,2);
%     v = zeros(2*n+1,nd);
%     Nt = size(x_t,1);
%     xf = fft(x_t);
%     v(1,:) = real(xf(1,:))/Nt;
%     for i=1:n
%         hi = h(i+1);
%         v(2*i,:) = real(xf(hi+1,:))/(Nt/2);
%         v(2*i+1,:) = -imag(xf(hi+1,:))/(Nt/2);
%     end
    
    zi = uint32(find(h==0));
    h([1 zi]) = h([zi 1]);
    n = uint32(length(h));
    nd = uint32(size(x_t,2));
    if ~isempty(zi)
        n = n-1;
        siz = 2*n+1;
    else
        siz = 2*n;
    end
    v = zeros(siz,nd, class(x_t));
    Nt = size(x_t,1);
    xf = fft(x_t);
    if ~isempty(zi)
        v(1,:) = real(xf(1,:))/Nt;
        zi = 1;
    else
        zi = 0;
    end
    for i=1:n
        hi = h(i+zi);
        v(2*i-1+zi,:) = real(xf(hi+1,:))/(Nt/2);
        v(2*i+zi,:) = -imag(xf(hi+1,:))/(Nt/2);
    end
end
