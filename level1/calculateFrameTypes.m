function [count_ols, count_esh, count_lss, count_lps]=calculateFrameTypes(y)

    %Split into frames(278)
    frames= buffer(y(:,1), 2048, 1024, 'nodelay');

    [row_length, col_length]=size(frames);
    
    count_ols=2;
    count_esh=0;
    count_lps=0;
    count_lss=0;
    
    frame_types(1)="OLS";
    
    for i=2:col_length-1
        frame_types(i) = SSC(frames(1:2048,i), frames(1:2048,i+1), frame_types(i-1));
        if(frame_types(i)=="OLS")
            count_ols=count_ols+1;
        elseif(frame_types(i)=="ESH")
            count_esh=count_esh+1;
        elseif(frame_types(i)=="LSS")
            count_lss=count_lss+1;
        else
            count_lps=count_lps+1;
        end
    end
    
    frame_types(col_length)=frame_types(col_length-1);
    
end