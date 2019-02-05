function x = iAACoder3(AACSeq3, fNameOut)
    %Returns decoded

    sequence_lentgth = length(AACSeq3);
    xi = zeros(sequence_lentgth*1024,2);
    for i = 1:sequence_lentgth-1
        % TODO Huffman decode
        streamL = AACSeq3(i).chl.stream;
        sfcL = AACSeq3(i).chl.sfc;
        streamR = AACSeq3(i).chr.stream;
        sfcR = AACSeq3(i).chr.sfc;
        frameFL = iAACquantizer(streamL, sfcL, AACSeq3(i).chl.G, AACSeq3(i).frameType);
        frameFR = iAACquantizer(streamR, sfcR, AACSeq3(i).chr.G, AACSeq3(i).frameType);
        if (AACSeq3(i).frameType == "ESH")
            frameF=zeros(128,8,2);
            frameF(:,:,1) = iTNS(frameFL , AACSeq3(i).frameType, AACSeq3(i).chl.TNScoeffs );
            frameF(:,:,2) = iTNS(frameFR , AACSeq3(i).frameType, AACSeq3(i).chr.TNScoeffs );
            frameT = iFilterbank(frameF, AACSeq3(i).frameType, AACSeq3(i).winType);
        else
            frameF=zeros(1024,2);
            frameF(:,1) = iTNS(frameFL , AACSeq3(i).frameType, AACSeq3(i).chl.TNScoeffs);
            frameF(:,2) = iTNS(frameFR , AACSeq3(i).frameType, AACSeq3(i).chr.TNScoeffs);
            frameT = iFilterbank([frameF(:,1),frameF(:,2)], AACSeq3(i).frameType, AACSeq3(i).winType);
        end
        xi((i-1)*1024 + (1:2048),:) = xi((i-1)*1024 + (1:2048),:) + frameT;
    end

    % Write audio sequence to a file using 48 KHz
    audiowrite(fNameOut,xi,48000);

    if(nargout==1)
        x = xi;
    end
end
