function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)

table=load('TableB219.mat');
if (frameType =="ESH")
    short_fft=table.B219b;
    %Calculate audibility threshold
    for j = 1:8
        for n = 1:42
            P(n) = sum(frameF(short_fft(n,2)+1:short_fft(n,3)+1,j).^2);
            T(n) = P(n) / SMR(n,j);
        end
        S
        
    end
else
    long_fft=table.B219a;
    %Calculate audibility threshold
    for n = 1:69
        P(n) = sum(frameF(long_fft(n,2)+1:long_fft(n,3)+1).^2);
        T(n) = P(n) / SMR(n);
    end
    
    
    
end















end