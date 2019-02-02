function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

table=load('TableB219.mat');
if (frameType =="ESH")
    %Calculate audacity thresholds
    short_fft=table.B219b;
    
else
    %Calculate audacity thresholds
    long_fft=table.B219a;
    for n = 1:69
        P(n) = sum(frameF(long_fft(n,2)+1:long_fft(n,3)+1));
        T(n) = P(n) / SMR(n);
    end
    
end















end