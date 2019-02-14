function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
    % AACoder3 creates an encoded array [AACSeq3] of audio sequences coded by the advanced audio coding standard
    table=load('TableB219.mat');
    windowType = "KBD";
    huffLUT = loadLUT();
    forceCodebook = 12;
    audio = audioread(fNameIn);
    % Buffer the audio array using 2048 sample windows with 50 % overlap
    framesLeft = buffer([zeros(1024,1);audio(:,1);zeros(1024,1)], 2048, 1024, 'nodelay');
    framesRight = buffer([zeros(1024,1);audio(:,2);zeros(1024,1)], 2048, 1024, 'nodelay');
    
    sequence_length = size(framesLeft,2);
    % Use 1 full window (2048 samples) at the start of sequence  so that the
    % function of psychoacoustic model can be called to predict the first window 
    framesLeft = [zeros(2048,1) [zeros(1024,1);framesLeft(1:1024,1)] framesLeft];
    framesRight = [zeros(2048,1) [zeros(1024,1);framesRight(1:1024,1)] framesRight];
    
    AACSeq3(sequence_length) = struct();
    % First element of sequence
    for i = 1:sequence_length
        j = i + 2;
        if i ~= 1 && i ~= sequence_length
            AACSeq3(i).frameType = SSC([framesLeft(:,j) framesRight(:,j)],...
        [framesLeft(:,j+1) framesRight(:,j+1)], AACSeq3(i-1).frameType);
        elseif i == 1
            AACSeq3(i).frameType = "OLS";
        else
            if AACSeq3(sequence_length-1).frameType == "LSS"
                AACSeq3(sequence_length).frameType = "ESH";
            elseif AACSeq3(sequence_length-1).frameType == "LPS"
                AACSeq3(sequence_length).frameType = "OLS";
            else
                AACSeq3(sequence_length).frameType = AACSeq3(sequence_length-1).frameType;
            end
        end
        AACSeq3(i).winType = windowType;
        frameF = filterbank([framesLeft(:,j) framesRight(:,j)], AACSeq3(i).frameType, AACSeq3(i).winType);
        if (AACSeq3(i).frameType == "ESH")
            [ frameFoutL, TNScoeffsL ] = TNS(frameF(:,:,1), AACSeq3(i).frameType);
            [ frameFoutR, TNScoeffsR ] = TNS(frameF(:,:,2), AACSeq3(i).frameType);
        else
            [ frameFoutL, TNScoeffsL ] = TNS(frameF(:,1), AACSeq3(i).frameType);
            [ frameFoutR, TNScoeffsR ] = TNS(frameF(:,2), AACSeq3(i).frameType);
        end
        AACSeq3(i).chl.TNScoeffs = TNScoeffsL;
        SMR_L = psycho(framesLeft(:,j), AACSeq3(i).frameType, framesLeft(:,j-1), framesLeft(:,j-2));
        AACSeq3(i).chl.T = audibilityThr(frameFoutL, AACSeq3(i).frameType, SMR_L,  table);
        

        
        AACSeq3(i).chr.TNScoeffs = TNScoeffsR;
        SMR_R = psycho(framesRight(:,j),  AACSeq3(i).frameType, framesRight(:,j-1), framesRight(:,j-2));
        AACSeq3(i).chr.T = audibilityThr(frameFoutR, AACSeq3(i).frameType, SMR_R,  table);
        [S_L, sfcL, G_L] = AACquantizer(frameFoutL, AACSeq3(i).frameType, SMR_L);
        [S_R, sfcR, G_R] = AACquantizer(frameFoutR,  AACSeq3(i).frameType, SMR_R);
        if (AACSeq3(i).frameType == "ESH")
            sfcL = reshape(sfcL,328,1);
            sfcR = reshape(sfcR,328,1);
        end
        AACSeq3(i).chl.G = G_L;
        [AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook] = encodeHuff(S_L, huffLUT);
        [AACSeq3(i).chl.sfc, ~] = encodeHuff(sfcL, huffLUT, forceCodebook);
        AACSeq3(i).chr.G = G_R;
        [AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook] = encodeHuff(S_R, huffLUT);
        [AACSeq3(i).chr.sfc, ~] = encodeHuff(sfcR, huffLUT, forceCodebook);
        
    end

    save(fnameAACoded,'AACSeq3');
    
end

function T = audibilityThr(frameF, frameType, SMR,  table)
if (frameType =="ESH")
    short_fft=table.B219b;
    T = zeros(42,8);
    for j = 1:8
        % Calculate audibility threshold for each frequency band of the
        % current subframe
        P = zeros(42,1);
        for n = 1:42
            P(n) = sum(frameF(short_fft(n,2)+1:short_fft(n,3)+1,j).^2);
        end
        T(:,j) = P ./ SMR(:,j);
    end
else
    long_fft=table.B219a;
    % Calculate audibility threshold for each frequency band
    P = zeros(69,1);
    for n = 1:69
        P(n) = sum(frameF(long_fft(n,2)+1:long_fft(n,3)+1).^2);
    end
    T = P ./ SMR;
end
end