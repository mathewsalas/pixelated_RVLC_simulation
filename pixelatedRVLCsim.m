clear all; clf;
##variables 
fc = 1;
fs = 100;
t = 1/fs:1/fs:1;
##transmitter
ndat = 100000; #total number of bits
hdat = ndat/2; #half the total number of bits
b = [0 1 0 1 randi([0 1], 1, hdat); 0 0 1 1 randi([0 1], 1, hdat)] * 2 - 1;  #nrz data
u = (b(1,:) + j*b(2,:))/sqrt(2); #complex waveform
x = real(vec(u) * e.^(2*pi*i*fc*t)); #modulation
x = vec(x');                     #line up symbols
##channel
#quantization
nshutter = 8;     #number of shutters
BER = [];
for nshutter = 1:40
nrange   = max(x); #range of signal
q = round(nshutter*(x+nrange)/(2*nrange)); #qunatize the signal
#noise
EbNodB = 6.6;            #Eb/No in decibles
EbNo   = 10^(EbNodB/10); #conversion from dB to fraction
Eb     = 1/2;            #energy per bit
No     = Eb/EbNo;        #Eb/(Eb/No) = NoEb/Eb = No
L      = length(x);      #each sample of the signal recieves noise
sigma = sqrt(No/2);      #deviation(square root of varriance)
n = sigma*randn(L,1);    #noise generation
##reciever
#y = reshape((q-nshutter/2)*(2/nshutter).+n, [fs, length(x)/fs])'; #seperate symbols
y = reshape(((2/nshutter)*(q) - 1).+n, [fs, length(x)/fs])'; #seperate symbols
a = y .* e.^(-2*pi*i*fc*t);             #demodulate
c = a(:,fs/2) + a(:,fs/4);              #recollect data
##Analysis
p = zeros(1,length(c));      #empty vector
p(find(real(c) > 0)) += 1;   #real determination boundary
p(find(imag(c) > 0)) += j;   #imaginary determination boudary
p = (p*2-(1+j))/sqrt(2);     #nrz encoding
nERR = length(find(p ~= u));  #number of bits in error
BER = [BER nERR/(ndat+4)] ;         #error rate
BER_t = (1/2)*erfc(sqrt(EbNo));
dV = ((2/nshutter)*(round(nshutter*(sqrt(Eb)+nrange)/(2*nrange)))-1)^2 - Eb;
BER_q = (1/2)*erfc(sqrt((Eb+dV)/No));
endfor
##figures
for i = 1:40
  if(rem(i, 8) == 0)
    pritnf("/n");
  endif
  printf("%f, ", BER(i));
endfor
printf("\n");