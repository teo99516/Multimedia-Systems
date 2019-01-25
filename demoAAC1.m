function SNR = demoAAC1(fNameIn, fNameOut)
% Demonstrates the operation of an AAC encoder level 1
%

AACSeq1 = AACoder1(fNameIn);

decodedAudio = iAACoder1(AACSeq1, fNameOut);

SNR = snr(audio,audio - decodedAudio);
end