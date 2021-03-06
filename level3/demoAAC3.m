function [SNR, bitrate, compression] = demoAAC3(fNameIn, fNameOut, fNameAACoded)
    % Demonstrates the operation of an AAC encoder using components of
    % level 3 and both previous levels

    audio = audioread(fNameIn);

    tic;
    AACSeq3 = AACoder3(fNameIn, fNameAACoded);
    fprintf('Time elapsed for encoding of AAC Sequence is %f seconds\n',toc);
    
    % Calculate the size of the encoded sequence in bits
    bitsSize = 0;
    for i=1:length(AACSeq3)
        if (AACSeq3(i).frameType == "ESH")
            bitsSize = bitsSize+2*4*8*4; %TNSCoeffs (4 bits each)
            bitsSize = bitsSize+2*8*64; %G (64 bits each)
            bitsSize = bitsSize+2*8*16; %Codebook (16 bits each)      
        else
            bitsSize = bitsSize+2*4*4; %TNSCoeffs
            bitsSize = bitsSize+2*64; %G
            bitsSize = bitsSize+2*16; %Codebook
        end
        bitsSize = bitsSize+length(AACSeq3(i).chl.stream);
        bitsSize = bitsSize+length(AACSeq3(i).chl.sfc);
        bitsSize = bitsSize+length(AACSeq3(i).chr.stream);
        bitsSize = bitsSize+length(AACSeq3(i).chr.sfc);
        bitsSize = bitsSize+3; %Window(1 bit) + FrameType(2 bits)
    end
    
    tic;
    decodedAudio = iAACoder3(AACSeq3, fNameOut);
    fprintf('Time elapsed for decoding of AAC Sequence is %f seconds\n',toc);
    % Calculate mean Bitrate
    timeSize = size(audio,1)/48000;
    bitrate = bitsSize / timeSize;
    % Calculate uncompressed audio Bitrate
    uncompressedBitrate = 48000*2*16; %48KHz * 2Channels * 16bits = 1536 Kbps
    compression = bitrate/uncompressedBitrate;
    error = audio - decodedAudio(1024+(1:size(audio,1)),:);
    % Indexing starts from sample 1025 in order to skip the first half of the
    % added first overlap window that is needed to restore the whole signal
    figure('Name','Level 3 Decoded Signal Error','NumberTitle','off');
    plot(error);
    title('Input-Decoded Signal Error')
    xlabel('Sample #');
    ylabel('Error');
    legend('Left Channel','Right Channel');
    SNR_L = snr(audio(:,1),error(:,1));
    SNR_R = snr(audio(:,2),error(:,2));
    SNR = snr(audio,error);
    fprintf('Left Channel SNR = %.4f dB\n',SNR_L);
    fprintf('Right Channel SNR = %.4f dB\n',SNR_R);
    fprintf('Uncompressed size = %.0d bits\n',uncompressedBitrate * timeSize);
    fprintf('Compressed size = %.0d bits\n',bitsSize);
    fprintf('Uncompressed Bitrate = %.6f Kbps\n',uncompressedBitrate/1000);
    fprintf('Compressed Mean Bitrate = %.6f Kbps\n',bitrate/1000);
    fprintf('Compression Ratio = %.3f %%\n',100*compression);
end