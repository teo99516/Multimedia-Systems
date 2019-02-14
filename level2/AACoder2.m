function AACSeq2 = AACoder2(fNameIn)

    windowType = "KBD";

    audio = audioread(fNameIn);
    
    framesLeft = buffer([zeros(1024,1);audio(:,1);zeros(1024,1)], 2048, 1024, 'nodelay');
    framesRight = buffer([zeros(1024,1);audio(:,2);zeros(1024,1)], 2048, 1024, 'nodelay');

    sequence_lentgth = size(framesLeft,2);

    AACSeq2(sequence_lentgth) = struct();
   
    for i = 1:sequence_lentgth
        if(i~=1 && i~=sequence_lentgth)
            AACSeq2(i).frameType = SSC([framesLeft(:,i) framesRight(:,i)],...
            [framesLeft(:,i+1) framesRight(:,i+1)], AACSeq2(i-1).frameType);
        elseif i==1
            AACSeq2(1).frameType = "OLS";
        else
            if AACSeq2(sequence_lentgth-1).frameType == "ESH"
                AACSeq2(sequence_lentgth).frameType = "LPS";
            elseif AACSeq2(sequence_lentgth-1).frameType == "LSS"
                AACSeq2(sequence_lentgth).frameType = "ESH";
            else
                AACSeq2(sequence_lentgth).frameType = "OLS";
            end
        end
        
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

end