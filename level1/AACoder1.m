function AACSeq1 = AACoder1(fNameIn)
% Returns the encoded AACSeq1 of a given input file as a struct
%
windowType = "KBD";
audio = audioread(fNameIn);

framesLeft = buffer([zeros(1024,1);audio(:,1);zeros(1024,1)], 2048, 1024, 'nodelay');
framesRight = buffer([zeros(1024,1);audio(:,2);zeros(1024,1)], 2048, 1024, 'nodelay');

sequence_lentgth = size(framesLeft,2);

AACSeq1(sequence_lentgth) = struct();

for i = 1:sequence_lentgth
    if(i~=1 && i~=sequence_lentgth)
        AACSeq1(i).frameType = SSC([framesLeft(:,i) framesRight(:,i)],...
        [framesLeft(:,i+1) framesRight(:,i+1)], AACSeq1(i-1).frameType);
    elseif i==1
        firstFrame = SSC([framesLeft(:,i) framesRight(:,i)],...
        [framesLeft(:,i) framesRight(:,i)], "OLS");
        if firstFrame == "LSS"
            AACSeq1(1).frameType = "ESH";
        else
            AACSeq1(1).frameType = "OLS";
        end
    else
        if AACSeq1(sequence_lentgth-1).frameType == "ESH"
            AACSeq1(sequence_lentgth).frameType = "LPS";
        elseif AACSeq1(sequence_lentgth-1).frameType == "LSS"
            AACSeq1(sequence_lentgth).frameType = "ESH";
        else
            AACSeq1(sequence_lentgth).frameType = "OLS";
        end
    end
    AACSeq1(i).winType = windowType;
    frameF = filterbank([framesLeft(:,i) framesRight(:,i)], AACSeq1(i).frameType, AACSeq1(i).winType);
    if (AACSeq1(i).frameType == "ESH")
        AACSeq1(i).chl.frameF = frameF(:,:,1);
        AACSeq1(i).chr.frameF = frameF(:,:,2);
    else
        AACSeq1(i).chl.frameF = frameF(:,1);
        AACSeq1(i).chr.frameF = frameF(:,2);
    end
end

end