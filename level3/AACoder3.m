function AACSeq3 = AACoder3(fNameIn, fnameAACoded))

    windowType = "KBD";

    audio = audioread(fNameIn);
    
    framesLeft = buffer(audio(:,1), 2048, 1024, 'nodelay');
    framesRight = buffer(audio(:,2), 2048, 1024, 'nodelay');
    
    sequence_lentgth = size(framesLeft,2);
    frames(:,:,1) = framesLeft;
    frames(:,:,2) = framesRight;
    
    AACSeq2(sequence_lentgth) = struct();
    frameF = zeros(1024,2);
    % First element of sequence
    AACSeq2(1).frameType = "OLS";
    AACSeq2(1).winType = windowType;
    frameF = filterbank([framesLeft(:,1) framesRight(:,1)],"OLS", windowType);
    [ frameFout, TNScoeffs ] = TNS(frameF(:,1), AACSeq2(1).frameType);
    AACSeq2(1).chl.TNScoeffs = TNScoeffs;
    AACSeq2(1).chl.frameF = frameFout;
    [ frameFout, TNScoeffs ] = TNS(frameF(:,2), AACSeq2(1).frameType);
    AACSeq2(1).chr.TNScoeffs = TNScoeffs;
    AACSeq2(1).chr.frameF = frameFout;
    for i = 2:sequence_lentgth - 1
        AACSeq2(i).frameType = SSC([framesLeft(:,i) framesRight(:,i)],...
        [framesLeft(:,i+1) framesRight(:,i+1)], AACSeq2(i-1).frameType);
        AACSeq2(i).winType = windowType;
        frameF = filterbank([framesLeft(:,i) framesRight(:,i)], AACSeq2(i).frameType, AACSeq2(i).winType);
        if (AACSeq2(i).frameType == "ESH")
            [ frameFout, TNScoeffs ] = TNS(frameF(:,:,1), AACSeq2(i).frameType);
            AACSeq2(i).chl.TNScoeffs = TNScoeffs;
            AACSeq2(i).chl.frameF = frameFout;
            [ frameFout, TNScoeffs ] = TNS(frameF(:,:,2), AACSeq2(i).frameType);
            AACSeq2(i).chr.TNScoeffs = TNScoeffs;
            AACSeq2(i).chr.frameF = frameFout;

        else
            [ frameFout, TNScoeffs ] = TNS(frameF(:,1), AACSeq2(i).frameType);
            AACSeq2(i).chl.TNScoeffs = TNScoeffs;
            AACSeq2(i).chl.frameF = frameFout;
            [ frameFout, TNScoeffs ] = TNS(frameF(:,2), AACSeq2(i).frameType);
            AACSeq2(i).chr.TNScoeffs = TNScoeffs;
            AACSeq2(i).chr.frameF = frameFout;
        end
        
    end
    % Last element of sequence
    AACSeq2(sequence_lentgth).frameType = "OLS";
    AACSeq2(sequence_lentgth).winType = windowType;
    frameF = filterbank([framesLeft(:,sequence_lentgth) framesRight(:,sequence_lentgth)],"OLS", windowType);
    [ frameFout, TNScoeffs ] = TNS(frameF(:,1), AACSeq2(sequence_lentgth).frameType);
    AACSeq2(sequence_lentgth).chl.TNScoeffs = TNScoeffs;
    AACSeq2(sequence_lentgth).chl.frameF = frameFout;
    [ frameFout, TNScoeffs ] = TNS(frameF(:,2), AACSeq2(sequence_lentgth).frameType);
    AACSeq2(sequence_lentgth).chr.TNScoeffs = TNScoeffs;
    AACSeq2(sequence_lentgth).chr.frameF = frameFout;



end