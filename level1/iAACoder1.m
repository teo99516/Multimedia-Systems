function x = iAACoder1(AACSeq1, fNameOut)
% Returns decoded
%

sequence_lentgth = length(AACSeq1);
xi = zeros((sequence_lentgth+1)*1024,2);
for i = 1:sequence_lentgth
    if (AACSeq1(i).frameType == "ESH")
        frameF(:,:,1) = AACSeq1(i).chl.frameF;
        frameF(:,:,2) = AACSeq1(i).chr.frameF;
        frameT = iFilterbank(frameF, AACSeq1(i).frameType, AACSeq1(i).winType);
    else
        frameT = iFilterbank([AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF], AACSeq1(i).frameType, AACSeq1(i).winType);
    end
    
    xi((i-1)*1024 + (1:2048),:) = xi((i-1)*1024 + (1:2048),:) + frameT;
end

% Write audio sequence to a file using 48 KHz
audiowrite(fNameOut,xi(1025:end-1024,:),48000);

if(nargout==1)
    x = xi;
end

end
