% Calculation of the RR interval and implementation of the KALMAN filter on 
% the test signal The time interval between the impulses in the test signal
% are computed and stored in the variable orginterval.
% A copy of the interval data is taken into the variable interval so that 
% the interval can be corrupted with the measurement noise for the noise 
% performance study of the filter. The variable intervalesti gives the
% estimated values from the filter. Q and V represent the process and 
% measurement noise covariance


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
