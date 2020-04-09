%%algorithm prototyping script

clc
clear
close all

%% settings
datafolder="sensor_tests";
data_files = dir(strcat(datafolder, "/sensors_recording*"));
cal_file = dir(strcat(datafolder, "/*calibration*"));

%% get calibration
[~,~,calib]=xlsread(strcat(cal_file.folder, "\", cal_file.name));


[cal_acc,cal_mag,cal_fs]=parseData(calib);
cal_acc=mean(cal_acc);
cal_acc=9.81*cal_acc/norm(cal_acc);
cal_mag=mean(cal_mag);
cal_mag=cal_mag/norm(cal_mag);

%%process

%step thru all
fdoms=[];
areas=[];
for i =1:length(data_files)
    [~,~,data]=xlsread(strcat(data_files(i).folder, "\", data_files(i).name));
    [acc,mag,fs]=parseData(data);
    %remove first and last two seconds- trial setup
    acc=acc(200:end-200,:);
    mag=mag(200:end-200,:);
    
    time_length=size(acc,1)/fs;
    t=0:1/fs:(1/fs)*(size(acc,1)-1);
    magNorms=sqrt(mag(:,1).^2+ mag(:,2).^2+ mag(:,3).^2);
    mag=mag./magNorms;
    
    %get translation to multiply 
    %https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    u = mag(:,:);                      % a and b must be column vectors. all orientations
    base = cal_mag;                      % original orientation
    R=[];
    I=[1,0,0;0,1,0;0,0,1];
    ssc=@(x)[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
    RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
    for j=1:length(u)
       a=u(j,:);
       R{j}=RU(a,base);%transposing takes inverse rotation       
    end
    %subtract rotated gravity
    for j=1:length(acc)
        acc(j,:)=acc(j,:)-(R{j}*cal_acc')';
    end
    
    time=t;
    %separate variables
    acc=acc-mean(acc);
    accX=acc(:,1);
    accY=acc(:,2);
    accZ=acc(:,3);
    
    %high/low pass to get rid of artifacts Y = bandpass(X,Fpass,Fs)
    hp=0.1;
    lp=20;
    accX = lowpass(accX,lp,fs);
    accY = lowpass(accY,lp,fs);
    accZ = lowpass(accZ,lp,fs);
    
    X=fft((accX+accY+accZ));
    n = length(accX);       % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power = abs(X).^2/n;    % power of the DFT
    
    new=floor(length(f)/2);
    f=f(1:new);
    power=power(1:new);
    plot(f,power)
    
    f_dom=f(power==max(power));
    
    
    figure()
    plot(time, accX)
    hold on
    plot(time, accY)
    plot(time, accZ)

    title('Acceleration Data Without Gravity')
    xlabel('time (s)')
    ylabel('Acceleration (m/s2)')
    legend({'x','y','z'})
    hold off

    %integrate twice for position. subtract average speed to discount drift
    %get speed
    vx=(1/fs)*cumtrapz(accX);
    vy=(1/fs)*cumtrapz(accY);
    vz=(1/fs)*cumtrapz(accZ);
    
%     hp=.25;
%     vx=highpass(vx,hp,fs);
%     vy=highpass(vy,hp,fs);
%     vz=highpass(vz,hp,fs);
%     %subtract drift speed
%     vx=detrend(vx);
%     vy=detrend(vy);
%     vz=detrend(vz);
    avg_vx=conv(vx,.01*ones(100,1),'same');
    avg_vy=conv(vy,.01*ones(100,1),'same');
    avg_vz=conv(vz,.01*ones(100,1),'same');
    vx=vx-avg_vx;
    vy=vy-avg_vy;
    vz=vz-avg_vz;
    %get position
    x=(1/fs)*cumtrapz(vx);
    y=(1/fs)*cumtrapz(vy);
    z=(1/fs)*cumtrapz(vz);
    
    hp=.25;
    x=highpass(x,hp,fs);
    y=highpass(y,hp,fs);
    z=highpass(z,hp,fs);
%     avg_x=conv(x,.01*ones(100,1),'same');
%     avg_y=conv(y,.01*ones(100,1),'same');
%     avg_z=conv(z,.01*ones(100,1),'same');
%     x=x-avg_x;
%     y=y-avg_y;
%     z=z-avg_z;
    

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
    
    figure()
    plot(time, x)
    hold on
    plot(time, y)
    plot(time, z)

    title('Position Data')
    xlabel('time (s)')
    ylabel('Displacement (m)')
    legend({'x','y','z'})
    hold off

    figure()
    colormap(copper)
    scatter3(x,y,z,3,time)
    i;
    
    
    
    %X=fft((x.*y.*z));
    X=fftn([x,y,z]);
    n = length(x);       % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power = abs(X).^2/n;    % power of the DFT
    
    new=floor(length(f)/2);
    f=f(1:new);
    power=power(1:new);
    plot(f,power)
    
    f_dom=f(power==max(power));
    
    dist=sum(  sqrt( diff(x).^2 + diff(y).^2 + diff(z).^2 )  ) - sqrt( (x(end)-x(1)).^2 + (y(end)-y(1)).^2 + (z(end)-z(1)).^2 );
    nCycles=time(end)*f_dom;
    avgCirc=dist/nCycles;
    a=sqrt(8/5)*avgCirc/pi;
    area = pi*a*a*.5;
    
    
    [cps, ints] = cycle_calculations(x,y,z,t);
    
    fprintf("file: %s \nf_dom: %d, cps: %d, area: %d\n\n", data_files(i).name, f_dom, cps, mean(ints) );
    fdoms=[fdoms, f_dom];
    areas=[areas, area];
    close all

end


function [cps, int] = cycle_calculations(x,y,z,time)
    smallest_dt=.25;
    fs=1/(time(2)-time(1));
    %calculate cycles per second in time period
    [vals,locs] = findpeaks(x);
    xCycles=sum( [true; diff(locs)/fs>smallest_dt]);
    xlocs=locs( [true; diff(locs)/fs>smallest_dt]);
    [vals,locs] = findpeaks(y);
    yCycles=sum(  [true; diff(locs)/fs>smallest_dt]);
    ylocs=locs(  [true; diff(locs)/fs>smallest_dt]);
    [vals,locs] = findpeaks(z);
    zCycles=sum( [true; diff(locs)/fs>smallest_dt]);
    zlocs=locs( [true; diff(locs)/fs>smallest_dt]);
%     [vals,locs] = findpeaks(x);
%     xCycles=sum((vals>0) & [true; diff(locs)/fs>smallest_dt]);
%     xlocs=locs((vals>0) & [true; diff(locs)/fs>smallest_dt]);
%     [vals,locs] = findpeaks(y);
%     yCycles=sum((vals>0) & [true; diff(locs)/fs>smallest_dt]);
%     ylocs=locs((vals>0) & [true; diff(locs)/fs>smallest_dt]);
%     [vals,locs] = findpeaks(z);
%     zCycles=sum((vals>0) & [true; diff(locs)/fs>smallest_dt]);
%     zlocs=locs((vals>0) & [true; diff(locs)/fs>smallest_dt]);
    
    cps=mean([xCycles,yCycles,zCycles])/time(end);
    %calculate pushing area in last n seconds
    %aproximate: numerical integration over radial frame. find center point
    %for each cycle and integrate change in r by change in theta
    for i=1:length(xlocs)-1
        xhere=x(xlocs(i):xlocs(i+1));
        yhere=y(xlocs(i):xlocs(i+1));
        zhere=z(xlocs(i):xlocs(i+1));
        here=[xhere,yhere,zhere];
        center=mean(here);
        radii=here-center;
        dists=sqrt(   (radii(:,1).^2)   +    (radii(:,2).^2)  +    (radii(:,3).^2));
        speed=diff(radii,1);
        swept = cumsum(sqrt(speed(:,1).^2+speed(:,2).^2+speed(:,3).^2));
        angles = 2*pi*(swept)'/max(swept);
        %angles=2*pi*(1:length(dists))/length(dists); %this is an approximation- figure this out. something to do with speed vector
        %angles=angles(1:length(dists));%
        int(i)=calc_polar_area(dists, angles);
    end
end
    
function int = calc_polar_area(dists, angles)
dists(length(dists)+1)=dists(1);
angles(length(angles)+1)=angles(1)+2*pi;
int=0;
for i=1:length(angles)-1
    dtheta=angles(i+1)-angles(i);
    h=(mean([dists(i),dists(i+1)]));
    int=int+.5*h*h*sin(dtheta);            
end
end

function [acc,mag,fs] = parseData (data)
    %dump the first 100 (~.2 seconds) and last second
    data=data(100:end,:);
    mags=data(strcmp(data(:,2),'MAG'),:);
    accs=data(strcmp(data(:,2),'ACC'),:);
    mags=cell2mat(mags(:,[1,3,4,5]));
    accs=cell2mat(accs(:,[1,3,4,5]));
    nMags=size(mags,1);
    nAccs=size(accs,1);
    modparam = floor(nAccs/nMags); 
    accs=accs(1:modparam:end,:);
    %correct sizes
    while (size(accs,1)>size(mags,1))
        accs=accs(1:end-1,:);
    end
    while (size(mags,1)>size(accs,1))
        mags=mags(1:end-1,:);
    end
    fs=size(accs,1)/((accs(end,1)-accs(1,1))/(10^9));
    acc=accs(:,[2,3,4]);
    mag=mags(:,[2,3,4]);

end
