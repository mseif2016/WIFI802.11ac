%function [xout] = scrambler(xin,initial_seq)
initial_seq=[1 0 0 0 0 0 0];
xin=[1 0 0 1 0 1 0 1 1 0 0 0 1 0];
      for i=1:127   %127 shift
        temp=xor (initial_seq(4),initial_seq(7));   %adding
        initial_seq=[temp initial_seq(1:end)];
      end
      initial_seq=initial_seq(1:127);
      
      data_length = length(xin);
      RepInt = floor(data_length/127);
      RepRem = mod(data_length,127);
      scram = [kron(ones(1,RepInt), initial_seq) initial_seq(1:RepRem)];
      xout=xor(xin,scram);
      
subplot(2,1,1);
stairs(xin);
subplot(2,1,2);
stairs(xout);
