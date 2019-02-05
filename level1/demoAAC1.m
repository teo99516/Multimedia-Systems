function SNR = demoAAC1(fNameIn, fNameOut)
% Demonstrates the operation of an AAC encoder level 1
%
audio = audioread(fNameIn);
tic;
AACSeq1 = AACoder1(fNameIn);
fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);
tic;

count_ols=0;
count_esh=0;
count_lps=0;
count_lss=0;

eshs = [];
for i=1:length(AACSeq1)
    if(AACSeq1(i).frameType=="OLS")
        count_ols=count_ols+1;
    elseif(AACSeq1(i).frameType=="ESH")
        count_esh=count_esh+1;
        eshs = [eshs i];
    elseif(AACSeq1(i).frameType=="LSS")
        count_lss=count_lss+1;
    elseif(AACSeq1(i).frameType=="LPS")
        count_lps=count_lps+1;
    end
end

decodedAudio = iAACoder1(AACSeq1, fNameOut);
fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);
error = audio(1025:size(decodedAudio,1)-1024,:) - decodedAudio(1025:size(decodedAudio,1)-1024,:);
%plot(error)
SNR_L = snr(audio(1025:size(decodedAudio,1)-1024,1),error(:,1));
SNR_R = snr(audio(1025:size(decodedAudio,1)-1024,2),error(:,2));
SNR = snr(audio(1025:size(decodedAudio,1)-1024,:),error);
fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
fprintf('Left Channel SNR = %.4f dB\n',SNR_R);
end