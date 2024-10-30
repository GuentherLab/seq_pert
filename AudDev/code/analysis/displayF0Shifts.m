% First define pHz and spHZ for each trial (e.g., spHz_6 =
% trialData(6).audapData.shiftedPitchHz). Then define x for each trial
% (e.g. x = [1:length(pHz_6)])

%JM - Improve this to me more automated

figure;
    title('D0 shifts');
    subplot(3,2,1);
    x = (1:length(pHz_6));
            plot(x,pHz_6,'b');
        hold on;
            plot(x,spHz_6,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 6');
    
    subplot(3,2,2,'align');
    x = (1:length(pHz_9));
            plot(x,pHz_9,'b');
        hold on;
            plot(x,spHz_9,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 9');
     
    subplot(3,2,3,'align');
    x = (1:length(pHz_11));
            plot(x,pHz_11,'b');
        hold on;
            plot(x,spHz_11,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 11');
        
     subplot(3,2,4,'align');
     x = (1:length(pHz_19));
            plot(x,pHz_19,'b');
        hold on;
            plot(x,spHz_19,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 19');
    
     subplot(3,2,5,'align');
     x = (1:length(pHz_21));
            plot(x,pHz_21,'b');
        hold on;
            plot(x,spHz_21,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 21');
        
figure;
    title('U0 Shifts');

    subplot(3,2,1);
    x = (1:length(pHz_2));
            plot(x,pHz_2,'b');
        hold on;
            plot(x,spHz_2,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 2');
    
    subplot(3,2,2,'align');
    x = (1:length(pHz_14));
            plot(x,pHz_14,'b');
        hold on;
            plot(x,spHz_14,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 14');
     
    subplot(3,2,3,'align');
    x = (1:length(pHz_26));
            plot(x,pHz_26,'b');
        hold on;
            plot(x,spHz_26,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 26');
        
     subplot(3,2,4,'align');
     x = (1:length(pHz_34));
            plot(x,pHz_34,'b');
        hold on;
            plot(x,spHz_34,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 34');
        
     subplot(3,2,5,'align');
     x = (1:length(pHz_40));
            plot(x,pHz_40,'b');
        hold on;
            plot(x,spHz_40,'r');
        xlim([150 700]);
        ylim([150 250]);
        title('Trial 34');
