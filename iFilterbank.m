function frameT = iFilterbank(frameF, frameType, winType)
% Returns iMDCT transformed input frameF 1024x2 as frameT
%

    if (frameType == "ESH")
        frames(:,:,1) = imdct4(frameF(:,:,1));
        frames(:,:,2) = imdct4(frameF(:,:,2));
    else
        frameT(:,1) = imdct4(frameF(:,1));
        frameT(:,2) = imdct4(frameF(:,2));
    end
    
    N = 2048;
    if( winType=="KBD")
        
        w=kaiser(1024+1,4*pi); 

        w_left_2048=zeros(1024,1);
        w_right_2048=zeros(1024,1);
        for n=1:(N/2)
            %Left an right KBN windows (w-left is the inverse of w_right)
           w_left_2048(n)=sqrt( sum(w(1:n) )/sum(w(1:N/2+1)) ) ;
           w_right_2048(1025-n)=sqrt( sum(w(1:n) )/sum(w(1:N/2+1)) ) ;
        end
        
        w=kaiser(128+1,6*pi);
        
        w_left_256=zeros(128,1);
        w_right_256=zeros(128,1);
        for n=1:(256/2)
            %Left an right KBN windows (w-left is the inverse of w_right)
           w_left_256(n)=sqrt( sum(w(1:n) )/sum(w(1:128+1)) ) ;
           w_right_256(129-n)=sqrt( sum(w(1:n) )/sum(w(1:128+1)) ) ;
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
        frameT(1:N/2,:)=frameT(1:N/2,:).*[w_left_2048 w_left_2048];
        frameT(N/2+1:N,:)=frameT(N/2+1:N,:).*[w_right_2048 w_right_2048];
    elseif (frameType=="LSS")
        frameT(1:1024,:)=frameT(1:1024,:).*[w_left_2048 w_left_2048];
        frameT(1473:1600,:)=frameT(1473:1600,:).*[w_right_256 w_right_256];
        frameT(1601:N,:)=0;
    elseif (frameType=="LPS")
        frameT(1:448,:)=0;
        frameT(449:576,:)=frameT(449:576,:).*[w_left_256 w_left_256];
        frameT(1025:2048,:)=frameT(1025:2048,:).*[w_right_2048 w_right_2048];
    else
        % Restore ESH sequence
        frameT=zeros(2048,2);
        % Aply small window weights to [frames] array
        frames = permute(frames,[1 3 2]);
        for i=1:8
            frames(1:128,:,i)=frames(1:128, :,i).*[w_left_256 w_left_256];
            frames(129:256,:,i)=frames(129:256, :,i).*[w_right_256 w_right_256];
            % Inverse buffer the [frames] array into [frameT] array
            frameT((449:704) + (i-1) * 128,:) = frameT((449:704) + (i-1) * 128,:) +...
                frames(:,:,i);
        end
        
        
    end

end

function y = imdct4(x)
% IMDCT4 Calculates the Modified Discrete Cosine Transform
%   y = imdct4(x)
%
%   x: input signal (can be either a column or frame per column)
%   y: IMDCT of x
%
%   Vectorize ! ! !

% ------- imdct4.m -----------------------------------------
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

% We need these for furmulas below
N     = flen;
M     = N/2;
twoN  = 2*N;
sqrtN = sqrt(twoN);

% We need this twice so keep it around
t = (0:(M-1)).';
w = diag(sparse(exp(-j*2*pi*(t+1/8)/twoN)));

% Pre-twiddle
t = (0:(M-1)).';
c = x(2*t+1,:) + j*x(N-1-2*t+1,:);
c = (0.5*w)*c;

% FFT for N/2 points only !!!
c = fft(c,M);

% Post-twiddle
c = ((8/sqrtN)*w)*c;

% Preallocate rotation matrix
rot = zeros(twoN,fnum);

% Sort
t = (0:(M-1)).';
rot(2*t+1,:)   = real(c(t+1,:));
rot(N+2*t+1,:) = imag(c(t+1,:)); 
t = (1:2:(twoN-1)).';
rot(t+1,:) = -rot(twoN-1-t+1,:);

% Shift
t = (0:(3*M-1)).';
y(t+1,:) =  rot(t+M+1,:);
t = (3*M:(twoN-1)).';
y(t+1,:) = -rot(t-3*M+1,:);
end