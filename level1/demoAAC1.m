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
error = audio - decodedAudio(1024+(1:size(audio,1)),:);
% Indexing starts from sample 1025 in order to skip the first half of the
% added first overlap window that is needed to restore the whole signal
figure('Name','Level 1 Decoded Signal Error','NumberTitle','off');
plot(error);
title('Input-Decoded Signal Error')
xlabel('Sample #');
ylabel('Error');
SNR_L = snr(audio(:,1),error(:,1));
SNR_R = snr(audio(:,2),error(:,2));
SNR = snr(audio, error);
fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
end