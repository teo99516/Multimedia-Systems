function AACSeq1 = AACoder1(fNameIn)
% Returns the encoded AACSeq1 of a given input file as a struct transformed input frameT 1024x2 as frameF
%
windowType = "KBD";
audio = audioread(fNameIn);

framesLeft = buffer(audio(:,1), 2048, 1024, 'nodelay');
framesRight = buffer(audio(:,2), 2048, 1024, 'nodelay');

[row_length, col_length] = size(framesLeft);



AACSeq1(col_length) = struct();
AACSeq1(1).frameType = "OLS";
for i = 2:row_length-1
    AACSeq1(i).frameType = SSC(framesLeft(i,:), framesLeft(i + 1,:), AACSeq1(i - 1).frameType);
    AACSeq1(i).winType = windowType;
    AACSeq1(i).chl.frameF = filterbank(framesLeft(i,:), AACSeq1(i).frameType, AACSeq1(i).winType);
    AACSeq1(i).chr.frameF = filterbank(framesRight(i,:), AACSeq1(i).frameType, AACSeq1(i).winType);
end