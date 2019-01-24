function frameType = SSC(frameT, nextFrameT, prevFrameType)
   
    %Apply a high pass filter
    %Apply difference equation of the filter
    NextFrameFiltered(1)= 0.7548*nextFrameT(1);  
    for i =2:length(nextFrameT) 
        NextFrameFiltered(i)=0.7458*(nextFrameT(i)-nextFrameT(i-1))+ 0.5095*NextFrameFiltered(i-1);
    end

    %Calculate the squared sum of the 128 values
    %8 values is calculated
    count=0;
    for i= (448+128+1):128:(length(NextFrameFiltered)-448)
        count=count+1;
       squaredSum(count)=0;
       for j=0:127
           squaredSum(count)=squaredSum(count)+NextFrameFiltered(i+j)^2;
       end
      
    end

    %Calculate the attack values
    nextframeType="NO ESS";
      for i=2:8
        previeusSum=0;
        for j=i-1:-1:1
            previeusSum=previeusSum+squaredSum(j);
        end
        mean_prev=previeusSum/(i-1);
        attackValue(i)=squaredSum(i)/mean_prev;
       
        if (attackValue(i)>10&&squaredSum(i)>0.001)
            nextframeType="ESH";
            break;
        end
    end

    %frameType:
        %OLS: ONLY_LONG_SEQUENCE
        %LSS: LONG_START_SEQUENCE
        %ESH: EIGHT_SHORT_SEQUENCE
        %LPS: LONG_STOP_SEQUENCE
    
    if (prevFrameType=="OLS")
        if(nextframeType=="ESH")
            frameType="LSS";
        else
            frameType="OLS";
        end
    elseif (prevFrameType=="ESH")
        if(nextframeType=="ESH")
            frameType="ESH";
        else
            frameType="LPS";
        end
    elseif (prevFrameType=="LSS")
        frameType="ESH";
    else
        frameType="OLS" ;
    end
        
end

