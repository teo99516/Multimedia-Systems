function SNR = demoAAC2(fNameIn, fNameOut)
    % Demonstrates the operation of an AAC encoder level 2

    audio = audioread(fNameIn);

    tic;
    AACSeq2 = AACoder2(fNameIn);
    fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);

    tic;
    decodedAudio = iAACoder2(AACSeq2, fNameOut);
    fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);

    error = audio(1025:size(decodedAudio,1)-1024,:) - decodedAudio(1025:size(decodedAudio,1)-1024,:);
    %plot(error)
    SNR_L = snr(audio(1025:size(decodedAudio,1)-1024,1),error(:,1));
    SNR_R = snr(audio(1025:size(decodedAudio,1)-1024,2),error(:,2));
    SNR = snr(audio(1025:size(decodedAudio,1)-1024,:),error);
    fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
    fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
    
end
