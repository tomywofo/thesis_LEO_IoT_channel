function k = throughput(Nbits, Ncode, BER)
k = Nbits/Ncode*(1-BER).^(Ncode);