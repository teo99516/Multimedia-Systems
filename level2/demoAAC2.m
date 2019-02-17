function SNR = demoAAC2(fNameIn, fNameOut)
% Demonstrates the operation of an AAC encoder using components of
% level 2 and level 1
% Inputs:
%     fNameIn   :   String that contains the name of the audio file that is
%                   about to be encoded according to AAC
%     fNameOut  :   String containing the name of the file that the Decoded
%                   audio signal will be saved into
% Output:
%     SNR       :   Signal to Noise Ratio between Input & Decoded
%                   audio signal in dB

audio = audioread(fNameIn);

tic;
AACSeq2 = AACoder2(fNameIn);
fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);

tic;
decodedAudio = iAACoder2(AACSeq2, fNameOut);
fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);

error = audio - decodedAudio(1024+(1:size(audio,1)),:);
% Indexing starts from sample 1025 in order to skip the first half of the
% added first overlap window that is needed to restore the whole signal
figure('Name','Level 2 Decoded Signal Error','NumberTitle','off');
plot(error);
title('Input-Decoded Signal Error')
xlabel('Sample #');
ylabel('Error');
SNR_L = snr(audio(:,1),error(:,1));
SNR_R = snr(audio(:,2),error(:,2));
SNR = snr(audio,error);
fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
end
