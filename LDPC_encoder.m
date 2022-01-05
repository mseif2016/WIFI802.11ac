function encoded_data = LDPC_encoder( data , A_MPDU_length , mstbc , MCS_index , BW , NSS)

% Kareem - Seif - Rafik 

MCS  = struct('Modulation' , {'BPSK' , 'QPSK' , 'QPSK' , '16-QAM', '16-QAM' , '64-QAM' , '64-QAM' , ...
                              '64-QAM' ,  '256-QAM' , '256-QAM' ,},...
              'Code_rate' , {0.5 , 0.5, 0.75, 0.5, 0.75, 2/3 ,3/4 , 5/6, 3/4, 5/6} ...
             , 'NBPSCS' , {1 , 2 , 2 , 4 ,4, 6, 6, 6, 8, 8});

R = MCS(MCS_index).Code_rate; 
Ncarriers = [52 , 108, 234 , 468];
NDBPS = Ncarriers(BW) * R * MCS(MCS_index).NBPSCS * NSS ;

%% Step one: OFDM Symbol number
Nsym = mstbc *  ceil((8 * A_MPDU_length +16)/(mstbc * NDBPS));
Npld = Nsym * NDBPS;

%% Step two : Code word size & Number of code words
NCBPS = NDBPS / R;
NTCB = NCBPS * Nsym;


if NTCB <= 648
    Ncw = 1;
    if NTCB >= Npld + 912 * (1-R)
        NLDPC = 1296;
    else
        NLDPC = 648;
    end
elseif NTCB <= 1296 & NTCB > 648
    Ncw = 1;
    if NTCB >= Npld + 1464 * (1 - R)
        NLDPC = 1944;
    else
        NLDPC = 1296;
    end
elseif NTCB <= 1944 & NTCB > 1296
    Ncw = 1;
    NLDPC = 1944;
elseif  NTCB  <= 2592 & NTCB > 1944 
    Ncw = 2;
    if NTCB >= Npld + 2916 * (1 - R)
    NLDPC = 1944;
    else 
    NLDPC = 1296;
    end
else
    NLDPC = 1944;
   Ncw = ceil(Npld/(1944*R));
end



%% Step Three : Shortening zero bits
Nshrt = max([0, Ncw * NLDPC * R - Npld]);
if Nshrt 
    Nspcw = floor(Nshrt/Ncw);
    shrt_array = zeros(1,Ncw);
    shrt_array(1:mod(Nshrt,Ncw)) = Nspcw +1;
    shrt_array(mod(Nshrt,Ncw)+1 : end) = Nspcw;
end

encoded = zeros(Ncw,NLDPC);
for i =1:Ncw
    temp = [data((i-1)*NLDPC*R + 1:i*NLDPC*R) zeros(1,shrt_array(i))];
    assert(length(temp)<=NLDPC,'Error in line 67, length of data&shtBits is bigger than NLDPC');
    encoded(i,1:length(temp)) = temp;
end
