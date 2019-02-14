function frameType = SSC(frameT, nextFrameT, prevFrameType)
   
    %Apply a high pass filter
    %Apply difference equation of the filter
    NextFrameFiltered = zeros(size(nextFrameT));
    NextFrameFiltered(1,:) = 0.7548*nextFrameT(1,:);  
    for i = 2:size(nextFrameT,1) 
        NextFrameFiltered(i,:)=0.7548*(nextFrameT(i,:)-nextFrameT(i-1,:))+ 0.5095*NextFrameFiltered(i-1,:);
    end

    %Calculate the squared sum of the 128 values
    %8 values is calculated
    count=0;
    squaredSum = zeros(2,8);
    for i= (448+128+1):128:(size(NextFrameFiltered,1)-448)
        count=count+1;
        squaredSum(:,count)=0;
        for j=0:127
            squaredSum(:,count)=squaredSum(:,count)+NextFrameFiltered(i+j,:).^2';
        end
      
    end

    %Calculate the attack values
    nextframeType(1:2)="NO ESS";
    for i=2:8
        previousSum = zeros(2,1);
        for j=i-1:-1:1
            previousSum(:)=previousSum(:)+squaredSum(:,j);
        end
        mean_prev(1:2)=previousSum(:)/(i-1);
        attackValue(i,:)=squaredSum(:,i)./(mean_prev');
       
        if (attackValue(i,1)>10&&squaredSum(1,i)>0.001)
            nextframeType(1) = "ESH";
            %break;
        end
        if (attackValue(i,2)>10&&squaredSum(2,i)>0.001)
            nextframeType(2) = "ESH";
            %break;
        end
    end

    %frameType: 
        %OLS: ONLY_LONG_SEQUENCE -> 1
        %LSS: LONG_START_SEQUENCE -> 2
        %ESH: EIGHT_SHORT_SEQUENCE -> 3
        %LPS: LONG_STOP_SEQUENCE -> 4
    channelFrameType(1:2) = "";
    for channel = 1:2
        if (prevFrameType=="OLS")
            if(nextframeType(channel)=="ESH")
                channelFrameType(channel)="LSS";
            else
                channelFrameType(channel)="OLS";
            end
        elseif (prevFrameType=="ESH")
            if(nextframeType(channel)=="ESH")
                channelFrameType(channel)="ESH";
            else
                channelFrameType(channel)="LPS";
            end
        elseif (prevFrameType=="LSS")
            channelFrameType(channel)="ESH";
        else
            channelFrameType(channel)="OLS" ;
        end
    end
    % Decide for current frameType based on each channel's seqType
    % according to Table 1
    Table1 = ["OLS" "OLS" "OLS";"OLS" "LSS" "LSS";"OLS" "ESH" "ESH";"OLS" "LPS" "LPS";...
              "LSS" "OLS" "LSS";"LSS" "LSS" "LSS";"LSS" "ESH" "ESH";"LSS" "LPS" "ESH";...
              "ESH" "OLS" "ESH";"ESH" "LSS" "ESH";"ESH" "ESH" "ESH";"ESH" "LPS" "ESH";...
              "LPS" "OLS" "LPS";"LPS" "LSS" "ESH";"LPS" "ESH" "ESH";"LPS" "LPS" "LPS";];
    for Case=1:16
        if all(Table1(Case,1:2) == channelFrameType)
            frameType = Table1(Case,3);
            break;
        end
    end
    
end




