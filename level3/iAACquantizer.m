function frameF = iAACquantizer(S, sfc, G, frameType)
table=load('TableB219.mat');
if (frameType =="ESH")
    % Initialize ESH
    short_fft=table.B219b;
    frameF = zeros(128,8);
    % Generate alpha indices
    alpha_sf = zeros(42,1);
    alpha_indices = ((0:127) >= short_fft(:,2)) & ((0:127) <= short_fft(:,3));
    alpha_indices = sum((1:42)*alpha_indices,1);
    
    for j = 1:8
        % Restore alpha values from DPCM
        alpha_sf(1,1) = G(j,1);
        for n = 1:41
            alpha_sf(n+1,1) = alpha_sf(n,1) + sfc(n,j);
        end
        % Restore MDCT values from Quantized vector S
        currentSubframe = S((j-1)*128+(1:128),1);
        frameF(:,j) = sign(currentSubframe) .* (abs(currentSubframe).^(4/3)) .* 2.^(alpha_sf(alpha_indices)/4);
    end
else
    % Initialize non ESH
    long_fft=table.B219a;
    % Generate alpha indices
    alpha_sf = zeros(69,1);
    alpha_indices = ((0:1023) >= long_fft(:,2)) & ((0:1023) <= long_fft(:,3));
    alpha_indices = sum((1:69)*alpha_indices,1);
    
    % Restore alpha values from DPCM
    alpha_sf(1,1) = G;
    for n = 1:68
        alpha_sf(n+1,1) = alpha_sf(n,1) + sfc(n,1);
    end
    % Restore MDCT values from Quantized vector S
    frameF = sign(S) .* (abs(S).^(4/3)) .* 2.^(alpha_sf(alpha_indices)/4);
end

end
