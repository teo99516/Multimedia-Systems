function frameFout = iTNS(frameFin, frameType, TNScoeffs)

    if (frameType=="ESH")
        frameFout=zeros(128,8);
        for i=1:8
            frameFout(:,i)=filter(1,[ 1; -TNScoeffs(:,i)],frameFin(:,i));
        end
    else
        %Apply inverse FIR
        frameFout=filter(1,[ 1; -TNScoeffs],frameFin);
    end

end