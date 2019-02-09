function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fNameAACoded)
    % Demonstrates the operation of an AAC encoder level 3

    audio = audioread(fNameIn);

    tic;
    AACSeq3 = AACoder3(fNameIn, fNameAACoded);
    fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);

    tic;
    decodedAudio = iAACoder3(AACSeq3, fNameOut);
    fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);

    error = audio(1025:size(decodedAudio,1)-1024,:) - decodedAudio(1025:size(decodedAudio,1)-1024,:);
    plot(error)
    SNR_L = snr(audio(1025:size(decodedAudio,1)-1024,1),error(:,1));
    SNR_R = snr(audio(1025:size(decodedAudio,1)-1024,2),error(:,2));
    SNR2 = mySNR(audio(1025:size(decodedAudio,1)-1024,:),error);
    SNR = snr(audio(1025:size(decodedAudio,1)-1024,:),error);
    fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
    fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
    fprintf('Alternative SNR = %.4f dB\n',SNR2);
end


function SNR = mySNR(x,e)
    SNR = sum(sum((x./e).^2,2));
    SNR = sqrt(SNR/(size(x,1)*size(x,2)));
end