%------------------ Alamouti [2 x nRx]----------------------%
%Description :
%Applying STBC with 2 Tx antenna using BPSK modulation and
%generalized number of receivers,but preferable not be so much
%to avoid receiver design complexity.
%In this code we see the performance of STBC [2 x 1] & [2 x 2]
%configuration

tic

clear all
clc

%Constants and Initializers :

N=1e6;                                                                     %Number of data
Eb_N0=0:2:20;                                                              %Eb/No in dB                                                                     %Number of receivers

%Preallocation :

Simulated_Ber=zeros(1,length(Eb_N0));

%Transmitter :

data=rand(1,N)>0.5;                                                        %Random data
Symbol=2*data-1;                                                           %BPSK modulation

%STBC :

Alamouti_scheme=1/sqrt(2)*kron(reshape(Symbol,2,N/2),ones(1,2)) ;


for nRx=1:2

for Eb_No_range= 1:length(Eb_N0)
    
    %Channel Parameters :
    
    channel_coefficient = 1/sqrt(2)*(randn(nRx,N) + 1i*randn(nRx,N));      %Rayleigh channel(generated by complex gaussian dist.)
    AWGN_channel = 1/sqrt(2)*(randn(nRx,N) + 1i*randn(nRx,N));             %White gaussian noise, with unity power
    
    %Preallocation(needed in the for loop) :
    
    received_signal = zeros(nRx,N);
    yMod = zeros(nRx*2,N);
    hMod = zeros(nRx*2,N);
    hEq= zeros(2*nRx,N);
    
    for index= 1:nRx
        
        hMod = kron(reshape(channel_coefficient(index,:),2,N/2),ones(1,2));
        temp = hMod;
        hMod(1,(2:2:end)) = conj(temp(2,(2:2:end)));
        hMod(2,(2:2:end)) = -conj(temp(1,(2:2:end)));
        
        % Received signal :
        
        received_signal(index,:) = sum(hMod.*Alamouti_scheme,1) +...
            10^(-Eb_N0(Eb_No_range)/20)*AWGN_channel(index,:);
        
        % Receiver :
        
        yMod((2*index-1:2*index),:) = kron(reshape(received_signal(index,:),2,N/2)...
            ,ones(1,2));
        
        % Forming the equalization matrix :
        
        hEq((2*index-1:2*index),:) = hMod;
        hEq(2*index-1,(1:2:end)) = conj(hEq(2*index-1,(1:2:end)));
        hEq(2*index,(2:2:end)) = conj(hEq(2*index,(2:2:end)));
        
    end
    
    % Equalization Process:
    
    hEqPower = sum(hEq.*conj(hEq),1);
    estimated_symbols = sum(hEq.*yMod,1)./hEqPower;                                     
    estimated_symbols(2:2:end) = conj(estimated_symbols(2:2:end));                                    
    
    % Detection process (Hard decision)
    
    estimated_data = real(estimated_symbols)>0;
    
    %BER :
    
    Simulated_Ber(Eb_No_range)=biterr(data,estimated_data)/N;              %relative frequency approach              
end

%Plot :

if (nRx==1)
semilogy(Eb_N0,Simulated_Ber,'rx-');
elseif (nRx==2) 
semilogy(Eb_N0,Simulated_Ber,'bo-');
end
hold on
end

xlabel('Eb/No(dB)');
ylabel('BER');
legend('2 x 1', '2 x 2');
title('Alamouti Algorithm with BPSK modulation');
toc