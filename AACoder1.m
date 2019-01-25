function AACSeq1 = AACoder1(fNameIn)
% Returns the encoded AACSeq1 of a given input file as a struct transformed input frameT 1024x2 as frameF
%
windowType = "KBD";

audio = audioread(fNameIn);

framesLeft = buffer(audio(:,1), 2048, 1024, 'nodelay');
framesRight = buffer(audio(:,2), 2048, 1024, 'nodelay');

sequence_lentgth = size(framesLeft,2);
frames(:,:,1) = framesLeft;
frames(:,:,2) = framesRight;

AACSeq1(sequence_lentgth) = struct();
frameF = zeros(1024,2);
% First element of sequence
AACSeq1(1).frameType = "OLS";
AACSeq1(1).winType = windowType;
AACSeq1(1).chl.frameF = filterbank(framesLeft(:,1),"OLS", windowType);
AACSeq1(1).chr.frameF = filterbank(framesRight(:,1),"OLS", windowType);
for i = 2:sequence_lentgth - 1
    AACSeq1(i).frameType = SSC(framesLeft(:,i), framesLeft(:,i + 1), AACSeq1(i - 1).frameType);
    AACSeq1(i).winType = windowType;
    frameF = filterbank(framesLeft(:,i), AACSeq1(i).frameType, AACSeq1(i).winType);
    if (AACSeq1(i).frameType == "ESH")
        AACSeq1(i).chl.frameF = frameF(:,:,1);
        AACSeq1(i).chr.frameF = frameF(:,:,1);
    else
        AACSeq1(i).chl.frameF = frameF(:,1);
        AACSeq1(i).chr.frameF = frameF(:,1);
    end
end
% Last element of sequence
AACSeq1(sequence_lentgth).frameType = "OLS";
AACSeq1(sequence_lentgth).winType = windowType;
AACSeq1(sequence_lentgth).chl.frameF = filterbank(framesLeft(:,1),"OLS", windowType);
AACSeq1(sequence_lentgth).chr.frameF = filterbank(framesRight(:,1),"OLS", windowType);

end