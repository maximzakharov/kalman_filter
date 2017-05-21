% defining the length of the data
nlen=3000;

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

% calculating the intervals
for index=1:numpul-1
    orginterval(index)=pulnum(index+1)-pulnum(index);
end
plen=length(orginterval); % length of the interval array
error=randn(1,plen);

% process and measurement noise covariance
Q=0.5;
V=0.1;
error=error*sqrt(V); % measurement noise that corrupts the RR interval

for index=1:numpul-1
    % interval(index)=orginterval(index)+error(index);
    interval(index)=orginterval(index);
end

kallen=index; % length of the KALMAN filter
% estimating the measured intervals using KALMAN filter
% preallocate memory for all arrays
intervalesti=zeros(1,kallen);
xaposteriori=zeros(1,kallen);
residual=zeros(1,kallen);
papriori=ones(1,kallen);
paposteriori=ones(1,kallen);
k=zeros(1,kallen);
% define the system
a=1; % State transition matrix
h=1; % Measurement matrix
% initial guesses for state and aposteriri covariance
xaposteriori_0=2;
paposteriori_0=1;
% predictor equations
xapriori(1)=a*xaposteriori_0;
intervalesti(1)=h*xapriori(1);
residual(1)=interval(1)-intervalesti(1);
papriori(1)=a*a*paposteriori_0+Q;
% corrector equations
k(1)=h*papriori(1)/(h*h*papriori(1)+V);
paposteriori(1)=papriori(1)*(1-h*k(1));
xaposteriori(1)=xapriori(1)+k(1)*residual(1);
% calculating the rest of the values

for j=2:kallen
    % predictor equations
    xapriori(j)=a*xaposteriori(j-1);
    intervalesti(j)=h*xapriori(j);
    residual(j)=interval(j)-intervalesti(j);
    papriori(j)=a*a*paposteriori(j-1)+Q;
    % corrector equations
    k(j)=h*papriori(j)/(h*h*papriori(j)+V);
    paposteriori(j)=papriori(j)*(1-h*k(j));
    xaposteriori(j)=xapriori(j)+k(j)*residual(j);
end
% end


% Reading the ECG files from the database
% Specify Data
PATH='./'; % path, where data are saved
HEADERFILE= '102.hea'; % header-file in text format
ATRFILE= '102.atr';
% attributes-file in binary format
DATAFILE='102.dat';
% data-file
SAMPLES2READ=15000; % number of samples to be read
% load header data
fprintf(1,'\\n$> WORKING ON %s ...\n', HEADERFILE);
signalh= fullfile(PATH, HEADERFILE);
fprintf(1,'\\n$> WORKING ON %s ...\n', signalh);
fid1=fopen(signalh,'r');
z=fgetl(fid1);


A= sscanf(z, '%*s %d %d %d',[1,3]);

nosig= A(1); % number of signals

fprintf(1,'\\n$> number of signals = %d\n', nosig);

sfreq=A(2); % sample rate of data
clear A;
for k=1:nosig
    z= fgetl(fid1);
    A= sscanf(z, '%*s %d %d %d %d %d',[1,5]);
    dformat(k)= A(1);
    % format; here only 212 is allowed
    gain(k)= A(2);
    % number of integers per mV
    bitres(k)= A(3);
    % bitresolution
    zerovalue(k)= A(4);
    % integer value of ECG zero point
    firstvalue(k)= A(5);
    % first integer value of signal (to test for errors)
end;
fclose(fid1);
clear A;

% load binary data
if dformat~= [212,212], error('this script does not apply binary formats different to 212.'); end;
signald= fullfile(PATH, DATAFILE);
% data in format 212
fid2=fopen(signald,'r');
A= fread(fid2, [3, SAMPLES2READ], 'uint8')';%matrix with 3 rows, each 8 bits long,=2*12bit
fclose(fid2);
M2H= bitshift(A(:,2), -4);
M1H= bitand(A(:,2), 15);
PRL=bitshift(bitand(A(:,2),8),9); % sign-bit
PRR=bitshift(bitand(A(:,2),128),5); % sign-bit
M( : , 1)= bitshift(M1H,8)+ A(:,1)-PRL;
M( : , 2)= bitshift(M2H,8)+ A(:,3)-PRR;
if M(1,:) ~= firstvalue, error('inconsistency in the first bit values'); end;
switch nosig
    case 2
        M( : , 1)= (M( : , 1)- zerovalue(1))/gain(1);
        M( : , 2)= (M( : , 2)- zerovalue(2))/gain(2);
        TIME=(0:(SAMPLES2READ-1))/sfreq;
    case 1
        M( : , 1)= (M( : , 1)- zerovalue(1));
        M( : , 2)= (M( : , 2)- zerovalue(1));
        M=M';
        M(1)=[];
        sM=size(M);
        sM=sM(2)+1;
        M(sM)=0;
        M=M';
        M=M/gain(1);
        TIME=(0:2*(SAMPLES2READ)-1)/sfreq;
        clear A M1H M2H PRR PRL;
        fprintf(1,'\\n$> LOADING DATA FINISHED \n');
        fprintf(1,'\\n$> ALL FINISHED \n');
        
        % QRS DETECTION
        x1=M(:,1);
        time=TIME;
        % initialization area
        ylp=zeros(1,length(x1));
        yhp=zeros(1,length(x1));
        yavg=zeros(1,length(x1));
        ymov=zeros(1,length(x1));
        yder=zeros(1,length(x1));
        yout=zeros(1,length(x1));
        thresholdf1=zeros(1,length(x1));
        thresholdf2=zeros(1,length(x1));
        thresholdi1=zeros(1,length(x1));
        thresholdi2=zeros(1,length(x1));
        spkf=zeros(1,500);
        npkf=zeros(1,500);
        spki=zeros(1,500);
        npki=zeros(1,500);
        period=max(time)/length(time);
        t=0;
        N=37;
        % this is the width of the integration window
        latency=30; % this is the delay for which a second signal can't be used as a QRS
        times=0;
        qmax=zeros(2,length(x1));
        wmax=zeros(2,length(x1));
        out=zeros(1,length(x1));
        fout=zeros(1,length(x1)+50);
        iout=zeros(1,length(x1));
        a=0;
        b=0;
        delay=0;
        k=2;
        %filter peak variable
        d=2;
        %integrator peak variable
        set=1; %these next 2 flags allow for finding the peaks
        let=1;
        spkf(1)=0;
        npkf(1)=0;
        spki(1)=0;
        npki(1)=0;
        RRsum=0;
        n=0;
        count=0;
        t=0;
        e=0;
        
        RRavg1=0;
        RRavg2=0;
        wcount=0;
        countp=0;
        intmax=0;
        filtmax=0;
        %this is the main loop of the program. the input is incrementally stepped
        %through the filter, derivative, integration, etc. the delays in
        %processing are taken into account using the if statements to make sure the
        %data is used at the correct time.
        for i=13:length(x1),
            %ylp is the output of the low pass filter
            ylp(i)=(x1(i)-2*x1(i-6)+x1(i-12)+2*ylp(i-1)-ylp(i-2));
            if i>32
                %yhp is the output of the bandpass filter
                yhp(i)=(yhp(i-1)-ylp(i)/32+ylp(i-16)-ylp(i-17)+ylp(i-32)/32);
            end
            if i>35
                %increments the search window
                filtmax=filtmax+1;
                if filtmax==50
                    qmax(1,k)=yhp(i);
                    qmax(2,k)=i;
                    for p=0:50
                        if yhp(i-p)>qmax(1,k)
                            %the max and location are stored in this variable. this
                            %max could be either signal or noise.
                            qmax(1,k)=yhp(i-p);
                            qmax(2,k)=i-p;
                        end
                    end
                    filtmax=0;
                    peakf=qmax;
                    spkf(1)=qmax(1,2);
                    %set=0;
                    if n<3
                        %if n is small the peaks are based on the previous peaks
                        if qmax(1,k)>qmax(1,k-1)
                            spkf(k)=.125*peakf(1,k)+.875*spkf(k-1);
                            npkf(k)=npkf(k-1);
                        else
                            npkf(k)=.125*peakf(1,k)+.875*npkf(k-1);
                            spkf(k)=spkf(k-1);
                        end
                    else
                        %once several peaks are determined, the sinal peak is
                        
                        %determined if it is within a certain range of a known peak
                        if (.85*qpk(n-1)<qmax(1,k))&(qmax(1,k)<1.15*qpk(n-1))
                            spkf(k)=.125*peakf(1,k)+.875*spkf(k-1);
                            npkf(k)=npkf(k-1);
                        else
                            npkf(k)=.125*peakf(1,k)+.875*npkf(k-1);
                            spkf(k)=spkf(k-1);
                        end
                    end
                    %determines the threshold value of this range
                    thresholdf1(k)=npkf(k)+.25*(spkf(k)-npkf(k));
                    thresholdf2(k)=.5*thresholdf1(k);
                    %determines if the signal of the specified range is greater or
                    %less than the threshold of that range and rewrites the output
                    %relative to the main i variable
                    for p=0:49
                        if yhp(i-p)>thresholdf1(k)%&(wmax(2,d-1)+latency)<i
                            fout(i-p)=1;
                        else
                            fout(i-p)=0;
                        end
                    end
                    %increment the peak variable
                    k=k+1;
                end
                if i==500+35
                    %this corresponds to 2seconds(the largest time needed for 1 beat
                    [a t]=max(qmax(1,:));
                    %finds maximum point after several max locations
                    n=1;
                    for m=1:length(qmax)
                        %qpk are the peaks that are within an allowable maximum,
                        %whereas qmax is all peaks which will include noise peaks
                        if (.85*a<qmax(1,m))&(qmax(1,m)<1.15*a)
                            qpk(1,n)=qmax(1,m);
                            qpk(2,n)=qmax(2,m);
                            %once the first RR is found, the other variable based
                            %on this can be started
                            if n>1
                                RR(n)=qpk(2,n)-qpk(2,n-1);
                                RRsum=RRsum+RR(n);
                                RRavg2=RRavg1;
                                RRlow=.92*RRavg2;
                                RRhigh=1.16*RRavg2;
                                RRmiss=1.67*RRavg2;
                            end
                            n=n+1;
                        end
                    end
                    
                end
                %once i is large enough, a the program can continulally search for
                %peaks and continually update the RR's
                if i>=500+35
                    %determines a range of acceptable peaks of qrs
                    if (.85*a<qmax(1,k))&(qmax(1,k)<1.15*a)
                        qpk(1,n)=qmax(1,k);
                        qpk(2,n)=qmax(2,k);
                        RR(n)=qpk(2,n)-qpk(2,n-1);
                        %this determines the averages if there are 8 or fewer RR
                        %intervals
                        if n<9
                            RRsum=RRsum+RR(n);
                            RRavg1=RRsum/(n-1);
                            RRavg2=RRavg1;
                        end
                        RRavg1=((RRavg1*.875)+.125*RR(n));
                        %this if statement determines RRbar and keeps track of how
                        %often RRbar is equal to RR.
                        if (RR(n)<RRhigh) & (RR(n)>RRlow)
                            RRbar(n)=RR(n);
                            RRavg2=.875*RRavg2+.125*RRbar(n);
                            count=count+1;
                            if count==8
                                RRavg2=RRavg1;
                            end
                        else
                            count=0;
                        end
                        RRlow=.92*RRavg2;
                        RRhigh=1.16*RRavg2;
                        RRmiss=1.67*RRavg2;
                        n=n+1;
                    end
                    %once several RR's are established, we can start looking for
                    %missed peaks
                    if n>3
                        if (i-qpk(2,n-1))>RRmiss
                            thresholdf1(k-1)=thresholdf2(k-1);
                            [peak t2]=max(peakf(1,qpk(2,n-2)+latency:i));
                            qpk(1,n-1)=peak;
                            qpk(2,n-1)=t2+qpk(2,n-2);
                            spkf(k-1)=0.25*peakf(1,k-1)+0.75*spkf(k-2);
                            RR(n-1)=qpk(2,n-1)-qpk(2,n-2);
                            if (RR(n-1)<RRhigh) & (RR(n-1)>RRlow)
                                RRbar(n-1)=RR(n-1);
                                RRavg2=.875*RRavg2+.125*RRbar(n-1);
                                count=count+1;
                                
                                if count==8
                                    RRavg2=RRavg1;
                                end
                            else
                                count=0;
                            end
                        end
                    end
                end
            end
            %this is the start of the integration process with the derivative and
            %the square outputs
            if i>36
                yder(i)=0.125*(2*yhp(i)+yhp(i-1)-yhp(i-3)-2*yhp(i-4));
                ysq(i)=yder(i)*yder(i);
            end
            %this is the start of the moving average integrator where is the size
            %of the window and yavg is the output that will be sent to the
            %threshold detector
            if i>72
                for j=0:N-1
                    ymov(i)=ymov(i)+ysq(i-j);
                end
                yavg(i)=ymov(i)/(N);
                %this section determines a maximum for a certain window(50samples)
                intmax=intmax+1;
                if intmax==50
                    wmax(1,d)=yavg(i);
                    wmax(2,d)=i;
                    for p=0:50
                        if yavg(i-p)>wmax(1,d)
                            %the max and location are stored in this variable. this
                            %max could be either signal or noise.
                            wmax(1,d)=yavg(i-p);
                            wmax(2,d)=i-p;
                        end
                    end
                    intmax=0;
                    %this initializer the first peak as a signal peak
                    spki(1)=wmax(1,2);
                    %let=0;
                    peaki=wmax;
                    thresholdi1(d)=npki(d-1)+.25*(spki(d-1)-npki(d-1));
                    thresholdi2(d)=.5*thresholdi1(d);
                    %average width of the integration window out from the threshold
                    %filter is equal to this delay
                    
                    delay=40;
                    %runs the average signal through the threshold at the same
                    %window interval that the threshold was determined
                    for p=0:49
                        if yavg(i-p)>thresholdi1(d)
                            iout(i-p)=1;
                        else
                            iout(i-p)=0;
                        end
                    end
                    %this for loop determines the output based on the delay. it is
                    %done here because we just calculated the iout(the most delayed
                    %signal).
                    for p=0:49
                        %we want to delay the fout signal by half the delay of the
                        %integration output
                        fdout(i-p)=fout(i-p-ceil(delay/2));
                        out(i-p)=fdout(i-p)*iout(i-p);
                    end
                    %this determines whether the peak is a signal or noise peak
                    if wmax(1,d)>wmax(1,d-1)
                        spki(d)=.125*peaki(1,d)+.875*spki(d-1);
                        npki(d)=npki(d-1);
                    else
                        npki(d)=.125*peaki(1,d)+.875*npki(d-1);
                        spki(d)=spki(d-1);
                    end
                    %the d variable is only incremented when there is a new peak
                    d=d+1;
                    %peak incrementer
                end
            end
            %once d and n are large enough. a missed peak in the integration area
            %can be determined
            if d>3 & n>3
                if wmax(2,d-2)+RRmiss<i
                    thresholdi1(d-1)=thresholdi2(d-1);
                    [peak t3]=max(peaki(1,wmax(2,d-2)+latency:i));
                    wmax(1,d-1)=peak;
                    wmax(2,d-1)=t3+wmax(2,d-2);
                    spki(d-1)=0.25*peaki(1,d-1)+0.75*spki(d-2);
                end
            end
        end
        
        % KALMAN algorithm implementation
        % calculating the RR Intervals from the QRS detector Output
        nlen=length(x1);
        numpul=0;
        index=4;
        while(index <= nlen)
            if (out(index-2)==0 && out(index)==1 && out(index-1)==0)
                numpul=numpul+1;
                pullocation(numpul)=index;
                index=index+160;
            else
                index=index+1;
            end
        end
        plen=length(pullocation);
        for index=1:plen-1;
            orginterval(index)=time(pullocation(index+1))-time(pullocation(index));
            interval(index)=time(pullocation(index+1))-time(pullocation(index));
        end
        kallen=plen-1;
        % estimating the measured intervals using kalman filter
        % preallocate memory for all arrays
        intervalesti=zeros(1,kallen);
        xaposteriori=zeros(1,kallen);
        residual=zeros(1,kallen);
        papriori=ones(1,kallen);
        paposteriori=ones(1,kallen);
        k=zeros(1,kallen);
        % define the system
        a=1; % State transition matrix
        h=1; % Measurement matrix
        Q=1; % Process noise covariance
        V=0.02; % Measurement noise covariance
        % initial guesses for state and aposteriri covariance
        xaposteriori_0=0.6;
        paposteriori_0=1;
        
        
        % predictor equations
        xapriori(1)=a*xaposteriori_0;
        intervalesti(1)=h*xapriori(1);
        residual(1)=interval(1)-intervalesti(1);
        papriori(1)=a*a*paposteriori_0+Q;
        % corrector equations
        k(1)=h*papriori(1)/(h*h*papriori(1)+V);
        paposteriori(1)=papriori(1)*(1-h*k(1));
        xaposteriori(1)=xapriori(1)+k(1)*residual(1);
        % calculating the rest of the values
        for j=2:kallen
            % predictor equations
            xapriori(j)=a*xaposteriori(j-1);
            intervalesti(j)=h*xapriori(j);
            residual(j)=interval(j)-intervalesti(j);
            papriori(j)=a*a*paposteriori(j-1)+Q;
            % corrector equations
            k(j)=h*papriori(j)/(h*h*papriori(j)+V);
            paposteriori(j)=papriori(j)*(1-h*k(j));
            xaposteriori(j)=xapriori(j)+k(j)*residual(j);
        end
end

