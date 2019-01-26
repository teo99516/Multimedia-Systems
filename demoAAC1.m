function SNR = demoAAC1(fNameIn, fNameOut)
% Demonstrates the operation of an AAC encoder level 1
%
audio = audioread(fNameIn);
tic;
AACSeq1 = AACoder1(fNameIn);
fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);
tic;
decodedAudio = iAACoder1(AACSeq1, fNameOut);
fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);
error = audio(1025:size(decodedAudio,1),:) - decodedAudio(1025:size(decodedAudio,1),:);
plot(error)
SNR_L = snr(audio(1025:size(decodedAudio,1)-1024,1),audio(1025:size(decodedAudio,1)-1024,1) - decodedAudio(1025:size(decodedAudio,1)-1024,1))
SNR_R = snr(audio(1025:size(decodedAudio,1)-1024,2),audio(1025:size(decodedAudio,1)-1024,2) - decodedAudio(1025:size(decodedAudio,1)-1024,2))
SNR = snr(audio(1025:size(decodedAudio,1)-1024,:),audio(1025:size(decodedAudio,1)-1024,:) - decodedAudio(1025:size(decodedAudio,1)-1024,:));
end