function pdelay=checkPertOn(ch1,ch2,timeVec,Fs,idxshift)

h1=figure;hold on;
plot(timeVec,ch1,'r')
plot(timeVec(1:end-idxshift+1),ch2(idxshift:end),'k');

h2=figure;
plot(timeVec(1:end-idxshift+1),ch1(1:end-idxshift+1)-ch2(561:end));

startPert = 1.01;
endRamp = 1.03;
shiftTime = idxshift/Fs;
rampTime = endRamp-startPert;
pdelay = startPert+shiftTime+rampTime


close(h1) 
close(h2)


