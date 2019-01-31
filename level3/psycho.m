function [r_T, phase_T] = psycho(frameT, frameType, frameTprev1, frameTprev2)

    % SMR 42X8 for ESH, else 69X1

    table=load('TableB219.mat');
  
    if(frameType=='ESH')
        short=zeros(42,42);
        % Calculate spreadingfunction for all band combinations of ESH frames
        for i =1:42
            for j=1:42
                short_fft=table.B219b;
                short(i,j) = spreadingfun(i, j, short_fft);
            end
        end
        for j=1:8
            N=length(frameT(:,j) );
            for n=1:256
                Sw=frameT(n,j)*( 0.5- 0.5*cos(pi*(n+5)/N) ) ;     
            end

            % Calculate FFT
            % Calculate abs and phase of every frequency(from 1 to 1024, the next are symmetrical) 
            fourier_trans=fft(Sw);
            r= abs(fourier_trans(1:128));
            phase= angle(fourier_trans(1:128));
        end
    else

        long= zeros(69,69);
        % Calculate spreadingfunction for all band combinations of non-ESH frames
        for i =1:69
            for j=1:69
                long_fft=table.B219a;
                long(i,j) = spreadingfun(i, j, long_fft);
            end
        end

        % Multiply with Hann windows
        N=length(frameT);
        for n=1:2048
            Sw_T(n)=frameT(n)*( 0.5- 0.5*cos(pi*(n+5)/N) ) ;
            Sw_prev1(n)=frameTprev1(n)*( 0.5- 0.5*cos(pi*(n+5)/N) ) ; 
            Sw_prev2(n)=frameTprev2(n)*( 0.5- 0.5*cos(pi*(n+5)/N) ) ; 
        end

        % Calculate FFT of this frame and the two before this
        % Calculate abs and phase of every frequency(from 1 to 1024, the next are symmetrical) 
        fourier_trans=fft(Sw_T);
        r_T= abs(fourier_trans(1:1024));
        phase_T= angle(fourier_trans(1:1024));

        fourier_trans=fft(Sw_prev1);
        r_prev1= abs(fourier_trans(1:1024));
        phase_prev1= angle(fourier_trans(1:1024));

        fourier_trans=fft(Sw_prev2);
        r_prev2= abs(fourier_trans(1:1024));
        phase_prev2= angle(fourier_trans(1:1024));
        
        % Calculate predicted abs and frequency
        r_T_pred=2*r_prev1-r_prev2;
        phase_pred=2*phase_prev1-phase_prev2;
        
        % Calculate predictability
        for i=1:2048
            temp1= r_T(n)*cos(phase_T(n))-r_T_pred(n)*cos(phase_pred());
            temp2= r_T(n)*sin(phase_T(n))- r_T_pred(n)*sin(phase_pred());
            
            predictability(n)= sqrt( (temp1)^2 + (temp2)^2 )/( r_T(n) + abs(r_T_pred) );
        end

        % Calculate energy and predictability for every band
        for i=1:69
            energy(i)= sum( r_T(long_fft(i,2)+1:long_fft(i,3)+1 )^2);
            predictability_2(i)= sum (predictability(long_fft(i,2)+1:long_fft(i,3)+1 )*r_T(long_fft(i,2)+1:long_fft(i,3)+1 )^2 );
        end

        % Combine energy and predictability 
        for n=1:69
            ecb(n)= sum(energy(1:69)*long(1:69, n));
            ct(n)= sum(predictability_2(1:69)*long(1:69, n));

            % Normalize predictability and energy 
            cb(n) = ct(n)/ecb(n);
            en(n) = ecb(n)/ sum(long(1:69, n) );

            % Calculate tonality index( values in (0,1) )
            tb(n) = -0.299 -0.43 *log(cb(n));

            % Noise Masking Tone( in dB)
            NMT(n) = 6;

            % Tone Masking Noise( in dB)
            TMN(n) = 18;

            %SNR
            SNR(n) = tb(n)*TMN(n) + (1 - tb(n))*NMT(n);

            % Convert from to Db to energy
            bc(n) = 10^(- SNR(n)/10);

            %Energy threshold
            nb(n) = en(n)*bc(n);





        end
       





        
        
        

    end



end

% Spreading function 
function x = spreadingfun(i, j, temp_table)
    
    if (i>=j)
        tmpx=3*(temp_table(j,5)- temp_table(i,5) );
    else
        tmpx = 1.5*(temp_table(j,5)- temp_table(i,5));
    end
    tmpz=8*min( (tmpx-0.5)^2-2*(tmpx-0.5),0);
    tmpy = 15.811389 + 7.5*(tmpx + 0.474) - 17.5*sqrt(1 + (tmpx + 0.474)^2 );
    if (tmpy < -100)
         x = 0;
    else
         x = 10^( (tmpz+tmpy)/10 );
    end

end