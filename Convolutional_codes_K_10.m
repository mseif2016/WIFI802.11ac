clc
clear all
%System Parameters:

k=10; %Input sequence length
uncodedWord=randi([0 1],1,k) %Input sequence

%Encoder:

L = 3;%length of shift register
polynomial=[7 3 5];%describing the connections of the encoder
trellis = poly2trellis(L,polynomial);
codeword=convenc(uncodedWord, trellis);



%QPSK Modulation:

codeword=reshape(codeword,2,length(codeword)/2);
modulated_symbols=(2*codeword(1,:))-1 + 1i*((2*codeword(2,:))-1);
modulated_symbols_modified=1/sqrt(2)*exp(1i*(pi/4))*modulated_symbols; 

%Performance at different SNR or Eb/No:
%Hint: BW=1/2T, then Eb/N0(dB)=SNR(dB)-3(dB)

 fc=5;
 T=1/fc;
 BW=1/2*T;
 for No=[2 4 8 12 16] %No={2,4,8,12,16}
    disp ('No = '); disp(No);
    received_sequence=modulated_symbols_modified...
        +No*BW*(randn(1,length(modulated_symbols_modified))+1i*randn(1,length(modulated_symbols_modified)));
    
    received_sequence_modified=sqrt(2)*exp(-1i*(pi/4))*received_sequence;
    
    Q_re = real(received_sequence_modified); %real
    Q_im = imag(received_sequence_modified); %imaginary
    
    
    
    %Bit detection:
    
    soft_bit_demod_signal=zeros(1,2*length(received_sequence));
    soft_bit_demod_signal(1:2:end)=Q_re;
    soft_bit_demod_signal(2:2:end)=Q_im;
    hard_bit_demodulation=soft_bit_demod_signal>0;
    
  
    %Viterbi decoder:
    
    %Hard decision decoding:
  
    traceback=2*L;
    recoveredWord=vitdec(hard_bit_demodulation,trellis,traceback,'trunc','hard')
      
    error = sum(abs(recoveredWord-uncodedWord));
    
    BER=error/k
 end
