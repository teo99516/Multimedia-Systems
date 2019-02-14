function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fNameAACoded)
    % Demonstrates the operation of an AAC encoder using components of
    % level 3 and both previous levels

    audio = audioread(fNameIn);

    tic;
    AACSeq3 = AACoder3(fNameIn, fNameAACoded);
    fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);
    
    bitsSum=0;
    for k=1:size(AACSeq3,1)
        if (strcmp(AACSeq3(k).frameType,'ESH'))
            for i=1:8
                bitsSum=bitsSum+size(AACSeq3(k).chl.stream{i},2);
                bitsSum=bitsSum+size(AACSeq3(k).chl.sfc{i},2);
                bitsSum=bitsSum+size(AACSeq3(k).chr.stream{i},2);
                bitsSum=bitsSum+size(AACSeq3(k).chr.sfc{i},2);
            end
            bitsSum=bitsSum+2*4*8*4; %TNSCoeffs( 4 bits each)
            bitsSum=bitsSum+2*8*64; %G (64 bits each)
            bitsSum=bitsSum+2*8*16; %Codebook (16 bits each)      
        else
            bitsSum=bitsSum+size(AACSeq3(k).chl.stream,2);
            bitsSum=bitsSum+size(AACSeq3(k).chl.sfc,2);
            bitsSum=bitsSum+size(AACSeq3(k).chr.stream,2);
            bitsSum=bitsSum+size(AACSeq3(k).chr.sfc,2);
            bitsSum=bitsSum+2*4*4; %TNSCoeffs
            bitsSum=bitsSum+2*64; %G 
            bitsSum=bitsSum+2*16; %Codebook
        end
        bitsSum=bitsSum+2*3*8; %Window + FrameType
    end

    tic;
    decodedAudio = iAACoder3(AACSeq3, fNameOut);
    fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);

    error = audio - decodedAudio(1024+(1:size(audio,1)),:);
    %plot(error)
    SNR_L = snr(audio(:,1),error(:,1));
    SNR_R = snr(audio(:,2),error(:,2));
    SNR = snr(audio,error);
    fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
    fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
    
end