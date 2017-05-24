% A train of impulses were generated. The first 1000 impulses have a time 
% interval of 10. The time interval for the next 1000 impulses was set to
% 30. The final 1000 impulses have a time interval of 10 again. The test
% signal was represented by the variable raw in the following code.


% defining the length of the data
nlen=3000;
% Loading the array initially with zeros
raw=zeros(1,nlen);
% generating the impulses with varying intervals
pulloc=1; % pulse location
i=1;
compre=10; % interval
numpul=0; % number of pulses
while(pulloc <= 1000)
    pulloc = pulloc + compre;
    width=0;
    while(width <= 1)
        if((width+pulloc)>nlen)
            break
        end
        raw(width+pulloc)=5;
        width=width+1;
    end
    if(pulloc <= nlen)
        pulnum(i)=pulloc;
        numpul=numpul+1;
    end
    i=i+1;
end
compre=30;
while(pulloc <= 2000)
    pulloc = pulloc + compre;
    width=0;
    while(width <= 1)
        if((width+pulloc)>nlen)
            break
        end
        raw(width+pulloc)=5;
        width=width+1;
    end
    if(pulloc <= nlen)
        pulnum(i)=pulloc;
        numpul=numpul+1;
    end
    i=i+1;
end
compre=10;
while(pulloc <= nlen)
    pulloc = pulloc + compre;
    width=0;
    while(width <= 1)
        if((width+pulloc)>nlen)
            break
        end
        raw(width+pulloc)=5;
        width=width+1;
    end
    if(pulloc <= nlen)
        pulnum(i)=pulloc;
        numpul=numpul+1;
    end
    i=i+1;
end


pulnum