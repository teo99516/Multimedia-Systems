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
        for n=1:(N)
            if(n<=N/2)
                w_left_2048(n)= sin( (pi/N)*(n+1/2)  );
            else
                w_right_2048(n-N/2)=sin( (pi/N)*(n+1/2)  );
            end
         end 
         for n=1:(256/2)
             if(n<=128)
                w_left_256(n)= sin( (pi/256)*(n+1/2)  );
            else
                w_right_256(n-128)=sin( (pi/256)*(n+1/2)  );
            end
        end
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
    frameF=frameT;
end