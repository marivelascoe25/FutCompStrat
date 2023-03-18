S = sparameters('passive.s2p');
freq = S.Frequencies;

tf_data = s2tf(S);
h = rational(freq,tf_data);

signalTime = [0,0.1,0.6,0.7,1.5]*1e-9;
signalValue = [0,5,5,0,0];
tper = 1.5e-9;

ts = 2e-11;
tsim = 0:ts:3*tper;
[tran,t] = pwlresp(h,signalTime,signalValue,tsim,tper);

vin = repmat(signalValue,1,3);
tin = [signalTime,signalTime+tper,signalTime+2*tper];
figure
plot(tin*1e9,vin,t*1e9,tran,'LineWidth',2)
axis([0 4.5 -2 6.2]);
xlabel('Time (ns)');
ylabel('Input Signal and Response (v)');
legend('Input','Resp');