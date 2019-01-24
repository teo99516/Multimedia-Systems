function frameF = filterbank(frameT, frameType, winType)
%winType: 
    %-KBN
    %-SIN
     N=length(frameT(:,1));
    %Kaiser windows
    %a=4 for N=2048 and a=6 for N=256
    if( winType=="KBN")
        
        w=kaiser(1024,4)  ; 

        w_left_2048=zeros(1024,1);
        w_right_2048=zeros(1024,1);
        for n=1:(N/2)
            %Left an right KBN windows (w-left is the inverse of w_right)
           w_left_2048(n)=sqrt( sum(w(1:n) )/sum(w(1:N/2)) ) ;
           w_right_2048(1025-n)=sqrt( sum(w(1:n) )/sum(w(1:N/2)) ) ;
        end
        
        w=kaiser(128,6);
        
        w_left_256=zeros(128,1);
        w_right_256=zeros(128,1);
        for n=1:(256/2)
            %Left an right KBN windows (w-left is the inverse of w_right)
           w_left_256(n)=sqrt( sum(w(1:n) )/sum(w(1:128)) ) ;
           w_right_256(129-n)=sqrt( sum(w(1:n) )/sum(w(1:128)) ) ;
        end
    else
        %Sinusoid windows
        w_left_2048=zeros(1024,1);
        w_right_2048=zeros(1024,1);
        for n=1:2048
            if(n<=1024)
                w_left_2048(n)= sin( (pi/N)*(n+1/2)  );
            else
                w_right_2048(n-N/2)=sin( (pi/N)*(n+1/2)  );
            end
        end 
         w_left_256=zeros(128,1);
        w_right_256=zeros(128,1);
         for n=1:(256/2)
             if(n<=128)
                w_left_256(n)= sin( (pi/256)*(n+1/2)  );
            else
                w_right_256(n-128)=sin( (pi/256)*(n+1/2)  );
            end
         end
        disp( length(w_left_2048) );
    end
    
    
    
    if (frameType=="OLS")
        frameT(1:N/2)=frameT(1:N/2).*w_left_2048;
        frameT(N/2+1:N)=frameT(N/2+1:N).*w_right_2048;
    elseif (frameType=="LSS")
        frameT(1:1024)=frameT(1:1024).*w_left_2048;
        frameT(1473:1600)=frameT(1473:1600).*w_right_256;
        frameT(1601:N)=0;
    elseif (frameType=="LPS")
        frameT(1:448)=0;
        frameT(449:576)=frameT(449:576).*w_left_256;
        frameT(1025:2048)=frameT(1025:2048).*w_right_2048;
    else
        %Keep the 1152 in the middle and make 8 subframes of 256 samples 
        frameT=frameT(449:1600);
        %Make a 265X8 array with the 8 subsamples
        frames= buffer(frameT, 256, 128, 'nodelay');
        for i=1:8
            frames(1:128,i)=frames(1:128, i).*w_left_256;
            frames(129:256,i)=frames(129:256, i).*w_right_256;
        end
    end
    % Calculate MDCT
    if (frameType == "ESH")
        N = 256;
        frameF = zeros(N/2,8);
%         n = 0:N-1;
%         k = 0:N/2-1;
%         cosineArgs = (2 * pi * ( n + (N/2 + 1)/2 ) / N)' * (k + 1/2);
        for i = 1:8
%             frameF(:,i) = 2 * sum(frames(:,i).*cos(cosineArgs));
            frameF(:,i) = mdct4(frames(:,i));
        end
    else
%         N = 2048;
%         n = 0:N-1;
%         k = 0:N/2-1;
%         frameF = 2 * sum(frameT .* cos((2 * pi * ( n + (N/2 + 1)/2 ) / N)' * k + 1/2))';
        frameF = mdct4(frameT);
    end
    
end

function y = mdct4(x)
% MDCT4 Calculates the Modified Discrete Cosine Transform
%   y = mdct4(x)
%
%   Use either a Sine or a Kaiser-Bessel Derived window (KBDWin)with 
%   50% overlap for perfect TDAC reconstruction.
%   Remember that MDCT coefs are symmetric: y(k)=-y(N-k-1) so the full
%   matrix (N) of coefs is: yf = [y;-flipud(y)];
%
%   x: input signal (can be either a column or frame per column)
%      length of x must be a integer multiple of 4 (each frame)     
%   y: MDCT of x (coefs are divided by sqrt(N))
%
%   Vectorize ! ! !

% ------- mdct4.m ------------------------------------------
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2002 by Columbia University.
% All rights reserved.
% ----------------------------------------------------------

[flen,fnum] = size(x);
% Make column if it's a single row
if (flen==1)
    x = x(:);
    flen = fnum;
    fnum = 1;
end
% Make sure length is multiple of 4
if (rem(flen,4)~=0)
    error('MDCT4 defined for lengths multiple of four.');
end

% We need these for furmulas below
N     = flen; % Length of window
M     = N/2;  % Number of coefficients
N4    = N/4;  % Simplify the way eqs look
sqrtN = sqrt(N);

% Preallocate rotation matrix
% It would be nice to be able to do it in-place but we cannot
% cause of the prerotation.
rot = zeros(flen,fnum);

% Shift
t = (0:(N4-1)).';
rot(t+1,:) = -x(t+3*N4+1,:);
t = (N4:(N-1)).';
rot(t+1,:) =  x(t-N4+1,:);
clear x;

% We need this twice so keep it around
t = (0:(N4-1)).';
w = diag(sparse(exp(-j*2*pi*(t+1/8)/N)));

% Pre-twiddle
t = (0:(N4-1)).';
c =   (rot(2*t+1,:)-rot(N-1-2*t+1,:))...
   -j*(rot(M+2*t+1,:)-rot(M-1-2*t+1,:));
% This is a really cool Matlab trick ;)
c = 0.5*w*c;
clear rot;

% FFT for N/4 points only !!!
c = fft(c,N4);

% Post-twiddle
c = (2/sqrtN)*w*c;

% Sort
t = (0:(N4-1)).';
y(2*t+1,:)     =  real(c(t+1,:));
y(M-1-2*t+1,:) = -imag(c(t+1,:));
end