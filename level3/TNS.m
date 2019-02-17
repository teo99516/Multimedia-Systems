function [ frameFout, TNScoeffs ] = TNS(frameFin, frameType)

%frameFin,frameFout: MDCT coefficiences: 128X8 for ESH, else 1024X1 

%TNS coefficiences For ESH 4X8, else 4X1

    table=load('TableB219.mat');
    long_fft=table.B219a;
    short_fft=table.B219b;
    %Calculate normalization coefficiences according to Table B219a
    if frameType=="ESH"
        %Calculate 42 values of band energies for every subframe
        normalization_coef=zeros(128,8);
        TNScoeffs = zeros(4,8);
        frameFout=zeros(128,8);
        for j=1:8
            for i=1:42
                %Calculate a new array for the squared values of frameFin, 
                %The selected values that are taken in account each time, are specified 
                %by column wlow(longfft(i,2)) and the following one wlow(long_fff(i+1,2) )
                %In each index add one for matlab compatibility
                band_energy = sum((frameFin(short_fft(i,2)+1:short_fft(i,3)+1, j) ).^2);
           
                
                %Callculate the normalization coefficiencies of the selected values
                %by taking the root square of band's energy
                normalization_coef( short_fft(i,2)+1:short_fft(i,3)+1, j )=sqrt(band_energy);
            end
            %Make values less steep  
            for i=127:-1:1
                normalization_coef(i,j)=( normalization_coef(i, j) + normalization_coef(i+1, j) )/2;
            end
            normalized_mdct=zeros(128,1);
            normalized_mdct(1)=frameFin(1,j)/normalization_coef(1,j);
            for i=2:128
                normalization_coef(i,j)=( normalization_coef(i, j) + normalization_coef(i-1, j) )/2;
                normalized_mdct(i)=frameFin(i,j)/normalization_coef(i,j);
            end 

            %Calculate the 4 optimized coeffieciencies a of the FIR prediction filter
            optimized_a_coef = lpc(normalized_mdct,4);
            
            optimized_a_coef=transpose(optimized_a_coef);

            %Quantize the FIR coefficiencies
            for i=2:5
                optimized_a_coef(i)=quanti( optimized_a_coef(i) );
            end
            
            %Apply FIR filter
            frameFout(:,j)=filter(optimized_a_coef,1,frameFin(:,j));
            flag_stable=isstable([1],optimized_a_coef);
            if ~flag_stable
                warning('Unstable TNS Coefficients filter');
            end
            
            %Return the positive values of a
            TNScoeffs(:,j)= -optimized_a_coef(2:5);    
            
        end
    else
        frameFout=zeros(1024,1);
        normalization_coef=zeros(1024,1);
        for i=1:69
            %Store band energy of the band(as before) 
            band_energy = sum((frameFin(long_fft(i,2)+1:long_fft(i,3)+1)).^2);
            
            %Callculate the normalization coefficiencies of the selected values
            %by taking the root square of band's energy
            normalization_coef( long_fft(i,2)+1:long_fft(i,3)+1 )=sqrt(band_energy);
        end 

        %Make values less steep and calculate the normalized MDCT coeffieciencies
        for i=1023:-1:1
            normalization_coef(i)=( normalization_coef(i) + normalization_coef(i+1) )/2;
        end
        normalized_mdct=zeros(1024,1);
        normalized_mdct(1)=frameFin(1)/normalization_coef(1);
        for i=2:1024
            normalization_coef(i)=( normalization_coef(i) + normalization_coef(i-1) )/2;
            normalized_mdct(i)=frameFin(i)/normalization_coef(i);
        end

        %Calculate the 4 optimized coeffieciencies a of the FIR prediction filter(a negative is calculated)
        optimized_a_coef = lpc(normalized_mdct,4);
        optimized_a_coef=transpose(optimized_a_coef);

        %Quantize the FIR coefficiencies
        for i=2:5
            optimized_a_coef(i)=quanti( optimized_a_coef(i) );
        end
        
        %Apply FIR filter               
        frameFout=filter(optimized_a_coef,1,frameFin);
        flag_stable=isstable(1,optimized_a_coef);
        if ~flag_stable
            warning('Unstable TNS Coefficients filter');
        end
        
        %Return the positive values of a
        TNScoeffs= -optimized_a_coef(2:5);
    end         

end

% A simplified quantizer that quantizes in the largest value of every partition
function [quantized_value]= quanti(in_value)
    quantized_value = max(min(floor((in_value+0.1)*10)/10-0.05,0.75),-0.75); %16 levels
    % Create the quantizer functon
    % in_value = -0.9:0.0001:0.9;
    % quantized_value = max(min(floor((in_value+0.1)*10)/10-0.05,0.75),-0.75);
    % figure
    % plot([in_value NaN zeros(1,length(in_value))],[zeros(1,length(in_value)) NaN in_value],'color','black')
    % hold on;
    % for i = -0.7:0.1:0.7
    %     plot(i*ones(1,2),[-0.9;0.9],':','MarkerSize',0.1,'color','black');
    % end
    % for i = -0.75:0.1:0.75
    %     plot([-0.9;0.9],i*ones(1,2),':','MarkerSize',0.1,'color','black');
    % end
    % plot(in_value,quantized_value,'LineWidth',1,'color','red')
    % title('4 Bit Quantizer - Dequantizer Figure')
    % xlabel('Input Value');
    % ylabel('Dequantized Value');
    % xticks(-0.7:0.1:0.7)
    % yticks(-0.75:0.1:0.75)
end