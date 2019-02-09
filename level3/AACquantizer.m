function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
MagicNumber = 0.4054;
S = zeros(1024,1);
table=load('TableB219.mat');
if (frameType =="ESH")
    % Initialize ESH
    short_fft=table.B219b;
    frameFS = zeros(size(frameF));
    frameFX_hat = frameFS;
    sfc = zeros(41,8);
    G = zeros(8,1);
    % Scale factors initial values
    alpha_sf = ones(42,8) * (16/3) * (ones(8,1) * log2((max(frameF).^(3/4))/8191));
    alpha_indices = ((0:127) >= short_fft(:,2)) & ((0:127) <= short_fft(:,3));
    alpha_indices = sum((1:42)*alpha_indices,1)';
    
    for j = 1:8
        % Calculate audibility threshold for each frequency band
        P = zeros(42,1);
        P_e = P;
        for n = 1:42
            P(n,1) = sum(frameF(short_fft(n,2)+1:short_fft(n,3)+1,j).^2);
        end
        T = P ./ SMR(:,j);
        
        % Initialize logical array inc. sfc(:,j) is already initialized with zeros()
        inc = 1 == 1;
        while(max(abs(sfc(:,j)))<=60 && any(inc))
            % Quantize MDCT values 
            frameFS(:,j) = sign(frameF(:,j)) .* floor((abs(frameF(:,j)).*2.^(-alpha_sf(alpha_indices)/4)).^(3/4) + MagicNumber);
            % Restore MDCT values from Quantized vector S
            frameFX_hat(:,j) = sign(frameFS(:,j)) .* (abs(frameFS(:,j)).^(4/3)) .* 2.^(alpha_sf(alpha_indices)/4);
            % Calculate the error Power of each frequency band
            for n = 1:42
                P_e(n) = sum((frameF(short_fft(n,2)+1:short_fft(n,3)+1,j) -... 
                              frameFX_hat(short_fft(n,2)+1:short_fft(n,3)+1,j)).^2);
            end
            % Increase the scale factors of bands that have a error Power
            % lower than the audibilty threshold
            inc = P_e < T;
            alpha_sf(inc,j) = alpha_sf(inc,j) + ones(sum(inc),1);
            % Update sfc DPCM values
            sfc(:,j) = alpha_sf(2:42,j) - alpha_sf(1:41,j);
        end
        % The while loop's statement has ceased to apply. The last increase
        % transaction must be undone
        alpha_sf(inc,j) = alpha_sf(inc,j) - ones(sum(inc),1);
        frameFS(:,j) = sign(frameF(:,j)) .* floor((abs(frameF(:,j)).*2.^(-alpha_sf(alpha_indices)/4)).^(3/4) + MagicNumber);
        sfc(:,j) = alpha_sf(2:42,j) - alpha_sf(1:41,j);
        G(j,1) = alpha_sf(1,j);
        % Store subframe S values in the full vector
        S((j-1)*128+(1:128),1) = frameFS(:,j);
    end
else
    % Initialize non ESH
    long_fft=table.B219a;
    % Scale factors initial values
    alpha_sf = ones(69,1) * (16/3) * log2((max(frameF)^(3/4))/8191);
    alpha_indices = ((0:1023) >= long_fft(:,2)) & ((0:1023) <= long_fft(:,3));
    alpha_indices = sum((1:69)*alpha_indices,1);
    
    % Calculate audibility threshold for each frequency band
    P = zeros(69,1);
    P_e = P;
    for n = 1:69
        P(n,1) = sum(frameF(long_fft(n,2)+1:long_fft(n,3)+1).^2);
    end
    T = P ./ SMR;
    
    % Initialize sfc for entry in while loop and logical array inc
    sfc = 0;
    inc = 1 == 1;
    while(max(abs(sfc))<=60 && any(inc))
        % Quantize MDCT values 
        S = sign(frameF) .* floor((abs(frameF).*2.^(-alpha_sf(alpha_indices)/4)).^(3/4) + MagicNumber);
        % Restore MDCT values from Quantized vector S
        X_hat = sign(S) .* (abs(S).^(4/3)) .* 2.^(alpha_sf(alpha_indices)/4);
        % Calculate the error Power of each frequency band
        for n = 1:69
            P_e(n) = sum((frameF(long_fft(n,2)+1:long_fft(n,3)+1) -... 
                          X_hat(long_fft(n,2)+1:long_fft(n,3)+1)).^2);
        end
        % Increase the scale factors of bands that have a error Power
        % lower than the audibilty threshold
        inc = P_e < T;
        alpha_sf(inc,1) = alpha_sf(inc,1) + ones(sum(inc),1);
        % Update sfc DPCM values
        sfc = alpha_sf(2:69) - alpha_sf(1:68);
    end
    % The while loop's statement has ceased to apply. The last increase
    % transaction must be undone
    alpha_sf(inc,1) = alpha_sf(inc,1) - ones(sum(inc),1);
    S = sign(frameF) .* floor((abs(frameF).*2.^(-alpha_sf(alpha_indices)/4)).^(3/4) + MagicNumber);
    sfc = alpha_sf(2:69) - alpha_sf(1:68);
    G = alpha_sf(1,1);
end

end
