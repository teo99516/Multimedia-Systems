function AACSeq3 = AACoder3(fNameIn, fnameAACoded)
    table=load('TableB219.mat');
    windowType = "KBD";

    audio = audioread(fNameIn);
    
    framesLeft = buffer(audio(:,1), 2048, 1024, 'nodelay');
    framesRight = buffer(audio(:,2), 2048, 1024, 'nodelay');
    
    sequence_lentgth = size(framesLeft,2);
    
    framesLeft = [zeros(2048,1) [zeros(1024,1);framesLeft(1:1024,1)] framesLeft];
    framesRight = [zeros(2048,1) [zeros(1024,1);framesRight(1:1024,1)] framesRight];
    
    AACSeq3(sequence_lentgth) = struct();
    % First element of sequence
    for i = 1:sequence_lentgth
        j = i + 2;
        if i ~= 1 && i ~= sequence_lentgth
            AACSeq3(i).frameType = SSC([framesLeft(:,j) framesRight(:,j)],...
        [framesLeft(:,j+1) framesRight(:,j+1)], AACSeq3(i-1).frameType);
        elseif i == 1
            AACSeq3(i).frameType = "OLS";
        else
            if AACSeq3(sequence_lentgth-1).frameType == "LSS"
                AACSeq3(sequence_lentgth).frameType = "ESH";
            elseif AACSeq3(sequence_lentgth-1).frameType == "LPS"
                AACSeq3(sequence_lentgth).frameType = "OLS";
            else
                AACSeq3(sequence_lentgth).frameType = AACSeq3(sequence_lentgth-1).frameType;
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
        [S_L, sfcL, G_L] = AACquantizer(frameFoutL, AACSeq3(i).frameType, SMR_L);
        AACSeq3(i).chl.G = G_L;
        AACSeq3(i).chl.T = audibilityThr(frameFoutL, AACSeq3(i).frameType, SMR_L,  table);
        % TODO Huffman encode
        AACSeq3(i).chl.sfc = sfcL;
        AACSeq3(i).chl.stream = S_L;
        
        AACSeq3(i).chr.TNScoeffs = TNScoeffsR;
        SMR_R = psycho(framesRight(:,j),  AACSeq3(i).frameType, framesRight(:,j-1), framesRight(:,j-2));
        [SR, sfcR, G_R] = AACquantizer(frameFoutR,  AACSeq3(i).frameType, SMR_R);
        AACSeq3(i).chr.G = G_R;
        AACSeq3(i).chr.T = audibilityThr(frameFoutR, AACSeq3(i).frameType, SMR_R,  table);
        % TODO Huffman encode
        AACSeq3(i).chr.sfc = sfcR;
        AACSeq3(i).chr.stream = SR;
        
    end

    save(fnameAACoded,'AACSeq3');
    
end

function T = audibilityThr(frameF, frameType, SMR,  table)
if (frameType =="ESH")
    short_fft=table.B219b;
    for j = 1:8
        % Calculate audibility threshold for each frequency band
        P = zeros(42,1);
        for n = 1:42
            P(n) = sum(frameF(short_fft(n,2)+1:short_fft(n,3)+1,j).^2);
        end
        T = P ./ SMR(:,j);
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