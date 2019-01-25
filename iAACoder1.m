function x = iAACoder1(AACSeq1, fNameOut)
% Returns MDCT transformed input frameT 1024x2 as frameF
%

sequence_lentgth = length(AACSeq1);
x = zeros(sequence_lentgth*1024,2);
for i = 1:sequence_lentgth - 2
    frameT = iFilterbank([AACSeq1(i).chl.frameF AACSeq1(i).chr.frameF], AACSeq1(i).frameType, AACSeq1(i).winType);
    x((i-1)*1024 + (1:2048),:) = x((i-1)*1024 + (1:2048),:) + frameT;
end

audiowrite(fNameOut,x,48000);

end
