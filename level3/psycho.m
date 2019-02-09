function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)
    persistent short long
    % SMR 42X8 for ESH, else 69X1

    table=load('TableB219.mat');
    if(frameType=='ESH')
        
        short_fft=table.B219b;
        % Calculate spreadingfunction for all band combinations of ESH frames
        if isempty(short)
            short=zeros(42,42);
            for i = 1:42
                for j = 1:42
                    short(i,j) = spreadingfun(i, j, short_fft);
                end
            end
        end
        frameT = buffer(frameT(449:1600),256,128,'nodelay');
        if(size(frameTprev1,2)>1)
            prevFrames = [frameTprev1(:,7:8) frameT(:,1:7)];
        else
            prevFrames = [frameTprev1(1217:1472) frameTprev1(1345:1600) frameT(:,1:7)];
        end
        
        for j=1:8
            N=length(frameT(:,j) );
            hann_256 = ( 0.5 - 0.5*cos(pi*((0:N-1)+0.5)/(N/2)) )';
            Sw_T = frameT(:,j) .* hann_256;   
            Sw_prev1 = prevFrames(:,j+1) .* hann_256; 
            Sw_prev2 = prevFrames(:,j) .* hann_256;
            
            % Calculate FFT
            % Calculate abs and phase of every frequency(from 1 to 128, the next are symmetrical) 
            fourier_trans=fft(Sw_T);
            r_T= abs(fourier_trans(1:N/2));
            phase_T= angle(fourier_trans(1:N/2));

            fourier_trans=fft(Sw_prev1);
            r_prev1= abs(fourier_trans(1:N/2));
            phase_prev1= angle(fourier_trans(1:N/2));

            fourier_trans=fft(Sw_prev2);
            r_prev2= abs(fourier_trans(1:N/2));
            phase_prev2= angle(fourier_trans(1:N/2));

            % Calculate predicted abs and frequency
            r_T_pred=2*r_prev1-r_prev2;
            phase_pred=2*phase_prev1-phase_prev2;

            % Calculate predictability
            for n=1:N/2
                temp1= r_T(n)*cos(phase_T(n))-r_T_pred(n)*cos(phase_pred(n));
                temp2= r_T(n)*sin(phase_T(n))- r_T_pred(n)*sin(phase_pred(n));

                predictability(n,1)= sqrt( (temp1)^2 + (temp2)^2 )/( r_T(n) + abs(r_T_pred(n)) );
            end
            % Calculate energy and predictability for every band
            for i=1:42
                energy(i)= sum( r_T(short_fft(i,2)+1:short_fft(i,3)+1 ).^2);
                predictability_2(i)= sum (predictability(short_fft(i,2)+1:short_fft(i,3)+1 ).*r_T(short_fft(i,2)+1:short_fft(i,3)+1 ).^2 );
            end
            % Combine energy and predictability 
            for n=1:42
                ecb(n)= sum(energy(1:42)*short(1:42, n));
                ct(n)= sum(predictability_2(1:42)*short(1:42, n));

                % Normalize predictability and energy 
                cb(n) = ct(n)/ecb(n);
                en(n) = ecb(n)/ sum(short(1:42, n) );

                % Calculate tonality index( values in (0,1) )
                tb(n) = max(min(-0.299 -0.43 *log(cb(n)),0.999999999),0.0000000001);

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

                %Pre-echo control
                qthr_hat(n) = eps()*(N/2)*10^(short_fft(n,6)/10);
                npart(n) = max( nb(n), qthr_hat(n));

                %Calculate SMR
                SMR(n,j) = energy(n)/npart(n);
            end
            
        end
        
    else
        long_fft=table.B219a;
        % Calculate spreadingfunction for all band combinations of non-ESH frames
        if isempty(long)
            long = zeros(69,69);
            for i =1:69
                for j=1:69
                    long(i,j) = spreadingfun(i, j, long_fft);
                end
            end
        end

        % Multiply with Hann windows
        N=length(frameT);
        hann_2048 = ( 0.5- 0.5*cos(pi*((0:N-1)+0.5)/(N/2)) )';
        
        Sw_T = frameT .* hann_2048;
        Sw_prev1 = frameTprev1 .* hann_2048; 
        Sw_prev2 = frameTprev2 .* hann_2048; 

        % Calculate FFT of this frame and the two before this
        % Calculate abs and phase of every frequency(from 1 to 1024, the next are symmetrical) 
        fourier_trans=fft(Sw_T);
        r_T= abs(fourier_trans(1:N/2));
        phase_T= angle(fourier_trans(1:N/2));

        fourier_trans=fft(Sw_prev1);
        r_prev1= abs(fourier_trans(1:N/2));
        phase_prev1= angle(fourier_trans(1:N/2));

        fourier_trans=fft(Sw_prev2);
        r_prev2= abs(fourier_trans(1:N/2));
        phase_prev2= angle(fourier_trans(1:N/2));
        
        % Copy symmetric half
%         r_T = [r_T r_T(end:-1:1)];
%         phase_T = [phase_T phase_T(end:-1:1)];
%         r_prev1 = [r_prev1 r_prev1(end:-1:1)];
%         phase_prev1 = [phase_prev1 phase_prev1(end:-1:1)];
%         r_prev2 = [r_prev2 r_prev2(end:-1:1)];
%         phase_prev2 = [phase_prev2 phase_prev2(end:-1:1)];
        
        % Calculate predicted abs and frequency
        r_T_pred=2*r_prev1-r_prev2;
        phase_pred=2*phase_prev1-phase_prev2;
        
        % Calculate predictability
        for n=1:N/2
            temp1= r_T(n)*cos(phase_T(n))-r_T_pred(n)*cos(phase_pred(n));
            temp2= r_T(n)*sin(phase_T(n))- r_T_pred(n)*sin(phase_pred(n));
            
            predictability(n,1)= sqrt( (temp1)^2 + (temp2)^2 )/( r_T(n) + abs(r_T_pred(n)) );
        end

        % Calculate energy and predictability for every band
        for i=1:69
            energy(i,1)= sum( r_T(long_fft(i,2)+1:long_fft(i,3)+1 ).^2);
            predictability_2(i,1)= sum (predictability(long_fft(i,2)+1:long_fft(i,3)+1 ).*r_T(long_fft(i,2)+1:long_fft(i,3)+1 ).^2 );
        end

        % Combine energy and predictability 
        for n=1:69
            ecb(n)= sum(energy(1:69).*long(1:69, n));
            ct(n)= sum(predictability_2(1:69).*long(1:69, n));

            % Normalize predictability and energy 
            cb(n) = ct(n)/ecb(n);
            en(n) = ecb(n)/ sum(long(1:69, n) );

            % Calculate tonality index( values in (0,1) )
            tb(n) =max(min(-0.299 -0.43 *log(cb(n)),0.999999999),0.0000000001);

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
            
            %Pre-echo control
            qthr_hat(n) = eps()*(N/2)*10^(long_fft(n,6)/10);
            npart(n) = max( nb(n), qthr_hat(n));
            
            %Calculate SMR
            SMR(n,1) = energy(n)/npart(n);
        end
        
    end
    
end

% Spreading function 
function x = spreadingfun(i, j, temp_table)
    
    if (temp_table(i,5)>=temp_table(j,5))
        tmpx = 3*(temp_table(j,5)- temp_table(i,5) );
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
