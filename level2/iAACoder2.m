function x = iAACoder2(AACSeq2, fNameOut)

    %Returns decoded

    sequence_lentgth = length(AACSeq2);
    xi = zeros(sequence_lentgth*1024,2);
    for i = 1:sequence_lentgth - 1
        if (AACSeq2(i).frameType == "ESH")
            frameF=zeros(128,8,2);
            frameF(:,:,1) = iTNS(AACSeq2(i).chl.frameF , AACSeq2(i).frameType, AACSeq2(i).chl.TNScoeffs );
            frameF(:,:,2) = iTNS(AACSeq2(i).chr.frameF , AACSeq2(i).frameType, AACSeq2(i).chr.TNScoeffs );
            frameT = iFilterbank(frameF, AACSeq2(i).frameType, AACSeq2(i).winType);
        else
            frameF=zeros(1024,2);
            frameF(:,1) = iTNS(AACSeq2(i).chl.frameF , AACSeq2(i).frameType, AACSeq2(i).chl.TNScoeffs);
            frameF(:,2) = iTNS(AACSeq2(i).chr.frameF , AACSeq2(i).frameType, AACSeq2(i).chr.TNScoeffs);
            frameT = iFilterbank([frameF(:,1),frameF(:,2)], AACSeq2(i).frameType, AACSeq2(i).winType);
           
        end

        xi((i-1)*1024 + (1:2048),:) = xi((i-1)*1024 + (1:2048),:) + frameT;
    end

    % Write audio sequence to a file using 48 KHz
    audiowrite(fNameOut,xi,48000);

    if(nargout==1)
        x = xi;
    end
end
