%%algorithm prototyping script

clc
clear
close all

%load in relevant data from trial
%trial begins with 5 seconds of still+level holding for calibration
%fs = 
[~,~,acc]=xlsread('accelerometer.csv');
[~,~,gyr]=xlsread('gyroscope.csv');
[~,~,mag]=xlsread('magneticfield.csv');
%len minus 1 bc header row is there

fs=floor((length(acc(:,1))-1)/(acc{end,1}/1000));
%go from ms to seconds
time=[acc{2:end,1};]/1000;


%visualize data
data={acc, gyr, mag};
nam={'acc', 'gyr', 'mag'};
for i=1:length(data)
   this=data{i}; 
   figure()
   title(nam{i})
   subplot(2,2,1)
   plot(time, [this{2:end,2}])
   title(['x ', nam{i}])

   subplot(2,2,2)
   plot(time, [this{2:end,3}])
   title(['y ', nam{i}])

   subplot(2,2,3)
   plot(time, [this{2:end,4}])
   title(['z ', nam{i}])
    
end

acc=cell2mat(acc(2:end,2:end));
gyr=cell2mat(gyr(2:end,2:end));
gyr=gyr-mean(gyr(1:20,:));
mag=cell2mat(mag(2:end,2:end));

%get relative orientation, normalize to starting values
xAngle=cumtrapz(gyr(:,1))/fs;
%xAngle=xAngle-xAngle(1);
yAngle=cumtrapz(gyr(:,2))/fs;
%yAngle=yAngle-yAngle(1);
zAngle=cumtrapz(gyr(:,3))/fs;
%zAngle=zAngle-zAngle(1);

%plot angles
if true
    figure()

    subplot(2,2,1)
    plot(time, xAngle)
    title(['x angle'])

    subplot(2,2,2)
    plot(time, yAngle)
    title(['y angle'])

    subplot(2,2,3)
    plot(time, zAngle)
    title(['z angle'])
end

%get rotation vector for each point
magnit=sqrt(mag(:,1).^2+ mag(:,2).^2+ mag(:,3).^2);
mag=mag./magnit;
mag=mag;
close all
%get translation to multiply 
%https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
u = mag(:,:);                      % a and b must be column vectors. all orientations
base = mag(1,:);                      % original orientation
R=[];
I=[1,0,0;0,1,0;0,0,1];
skewsym=@(x)[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
for i=2:length(u)
   a=u(i,:);
   v= cross(a,base);
   vx=skewsym(v);
   R{i}=I+vx+vx.*vx*(1-dot(a,base))/norm(v)^2; 
end

%rotate acceleration back to starting frame
for i=2:length(acc)
   %acc(i,:)=R{i}*acc(i,:)';
end
%subtract gravity
acc=acc-(mean(acc(1:100,:)));

%separate variables
accX=acc(:,1);
accY=acc(:,2);
accZ=acc(:,3);
%high/low pass to get rid of artifacts Y = bandpass(X,Fpass,Fs)
hp=0.3;
lp=10;
accX = bandpass(accX,[hp,lp],fs);
accY = bandpass(accY,[hp,lp],fs);
accZ = bandpass(accZ,[hp,lp],fs);


X=fft(accX);
n = length(accX);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(X).^2/n;    % power of the DFT
plot(f,power)
xlabel('Frequency')
ylabel('Power')



if false
    figure()
    subplot(2,2,1)
    plot(time, accX)
    title(['x acc corrected'])
    subplot(2,2,2)
    plot(time, accY)
    title(['y acc corrected'])
    subplot(2,2,3)
    plot(time, accZ)
    title(['z acc corrected'])
end

%integrate twice for position. subtract average speed to discount drift
%get speed
vx=(1/fs)*cumtrapz(accX);
vy=(1/fs)*cumtrapz(accY);
vz=(1/fs)*cumtrapz(accZ);
%subtract drift speed
vx=vx-mean(vx);
vy=vy-mean(vy);
vz=vz-mean(vz);
%get position
x=(1/fs)*cumtrapz(vx);
y=(1/fs)*cumtrapz(vy);
z=(1/fs)*cumtrapz(vz);
% hp=.01;
% x=highpass(x,hp,fs);
% y=highpass(y,hp,fs);
% z=highpass(z,hp,fs);

figure()

subplot(2,2,1)
plot(time, vx)
title(['x v'])
subplot(2,2,2)
plot(time, vy)
title(['y v'])
subplot(2,2,3)
plot(time, vz)
title(['z v'])


x=x(100:end);
y=y(100:end);
z=z(100:end);
time=time(100:end);
time=time-time(1);

figure()
subplot(2,2,1)
plot(time, x)
title(['x posn'])
subplot(2,2,2)
plot(time, y)
title(['y posn'])
subplot(2,2,3)
plot(time, z)
title(['z posn'])

figure()
colormap(copper)
scatter3(x,y,z,3,time)
%leave it here- good enough


%% goals- 
    %calculate cycles per second in time period
    [vals,locs] = findpeaks(x);
    xCycles=sum((vals>0) & [true; diff(locs)/20>.33]);
    xlocs=locs((vals>0) & [true; diff(locs)/20>.33]);
    [vals,locs] = findpeaks(y);
    yCycles=sum((vals>0) & [true; diff(locs)/20>.33]);
    ylocs=locs((vals>0) & [true; diff(locs)/20>.33]);
    [vals,locs] = findpeaks(z);
    zCycles=sum((vals>0) & [true; diff(locs)/20>.33]);
    zlocs=locs((vals>0) & [true; diff(locs)/20>.33]);
    
    cps=mean([xCycles,yCycles,zCycles])/time(end);
    %calculate pushing area in last 5 seconds
    %aproximate: numerical integration over radial frame. find center point
    %for each cycle and integrate change in r by change in theta







