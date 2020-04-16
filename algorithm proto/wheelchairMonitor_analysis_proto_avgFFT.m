%%algorithm prototyping script

clc
clear
close all

%% frequency data file naming convention: 
    %first number means area = [.15m^2, .20m^2]
    %second number means freq = [2Hz, 1.5Hz, 1.0Hz, 0.5Hz]
    %third numer means which replicate with these settings
    
%% settings
% %current best combination - method1 r^2 of .6
% detrend_v=false;
% 
% remove_gravity=true; %keep this
% remove_mean_acc=false;
% lowpass_acc = true;
% highpass_acc = true;
% 
% highpass_dist = false;
% highpass_v=false;
% remove_avg_posn = true;
% 
% hp_v = .20;
% hp_dist = .20;
% lp_acc=3;
% hp_acc=.25;
% 
% remove_avg_speed_stag=false;
% remove_mean_v=false;
% 
% %length of time to select from full series
% desired_dT = 5;

%tying this combination
%tying this combination

detrend_v=true;

remove_gravity=true; %keep this
remove_mean_acc=false;
lowpass_acc = true;
highpass_acc = false;

highpass_dist = false;
highpass_v=false;
remove_avg_posn = false;

hp_v = .20;
hp_dist = .20;
lp_acc=3;
hp_acc=.25;

remove_avg_speed_stag=false;
remove_mean_v=true;

%length of time to select from full series


%% select the data folder you are processing
%both of these contain their own calibration file
%datafolder="sensor_tests_area";
datafolder="sensor_tests_frequency";
%datafolder="sensor_tests_activities";

%% read file names
data_files = dir(strcat(datafolder, "/sensors_recording*"));
cal_file = dir(strcat(datafolder, "/*calibration*"));

%% get calibration file
[~,~,calib]=xlsread(strcat(cal_file.folder, "\", cal_file.name));
%parse calibration file
[cal_acc,cal_mag,cal_fs]=parseData(calib);
cal_acc=mean(cal_acc);
%set proper magnitude but use direction
cal_acc=9.81*cal_acc/norm(cal_acc);
cal_mag=mean(cal_mag);
cal_mag=cal_mag/norm(cal_mag);

%%process

%step thru all files
fdoms=[]; %dominant frequencies
areas1=[]; %swept areas
areas2=[]; %swept areas
areas3=[]; %swept areas
areas4=[]; %swept areas
for i =1:length(data_files)
    %read in data file
    [~,~,data]=xlsread(strcat(data_files(i).folder, "\", data_files(i).name));
    [acc,mag,fs]=parseData(data);
    desired_dT = 40*(256/fs);
    
    %TODO quaternion stuff.
%     qe = ecompass(acc, mag);
%     temp = euler(qe(1,:), 'ZYX', 'frame');
%         
%         yaw_c = temp(1);
%         pitch_c = temp(2);
%         roll_c = temp(3);
%     for frame=1:1:1000
%     
%         %Convert Quaternion to Roll Pitch Yaw
%         temp = euler(qe(frame,:), 'ZYX', 'frame');
%         
%         yaw = temp(1);
%         pitch = temp(2);
%         roll = temp(3);
%         %Visualize Data On Sphere Or any Other Objects
%         [x,y,z] = sphere;h = surf(x,y,z);axis('square'); 
%         title('AOL TEAM')
%         xlabel('x'); ylabel('y'); zlabel('z');
%         %Rotate Object
%         rotate(h,[1,0,0],(roll-roll_c)*180/pi)
%         rotate(h,[0,1,0],(pitch-pitch_c)*180/pi)
%         rotate(h,[0,0,1],(yaw-yaw_c)*180/pi+90)
%         view(0,0);
%         drawnow
%     end
    
    
    %remove first and last two seconds of data to ignore trial setup and stop time
    acc=acc(200:end-200,:);
    mag=mag(200:end-200,:);
    
    %pick random desired_dT seconds
    desiredL = ceil(fs*desired_dT);
    if size(acc,1)>desiredL
        start=randi(size(acc,1)-desiredL);
        acc=acc(start:start+desiredL,:);
        mag=mag(start:start+desiredL,:);
    end
    
    %make time vector
    t_length=size(acc,1)/fs;
    t=0:1/fs:(1/fs)*(size(acc,1)-1);
    %normalize magnetic fields to unit vectors. needs to be unit vector for gravity removal
    %part
    magNorms=sqrt(mag(:,1).^2+ mag(:,2).^2+ mag(:,3).^2);
    mag=mag./magNorms;
    
    %get translation to rotate gravity into new reference orientation for
    %each time point
    %https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    u = mag(:,:);                      % input
    base = cal_mag;                      % original orientation
    R=[];
    I=[1,0,0;0,1,0;0,0,1];
    ssc=@(x)[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ]; %skewsymmetric
    RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);  %rotation calc function
    for j=1:length(u)
       a=u(j,:); %get mag for this time point
       R{j}=RU(a,base); %transposing takes inverse rotation       
    end
    %subtract rotated gravity
    for j=1:length(acc)
        if remove_gravity
            acc(j,:)=acc(j,:)-(R{j}*cal_acc')';
            %disp(norm(R{j}*cal_acc'))
        end
    end  
    
    
    
    %separate variables
    if remove_mean_acc
        acc=acc-mean(acc);
    end
    
    accX=acc(:,1);     accY=acc(:,2);    accZ=acc(:,3);
    
      
    %do fft frequency analysis
    %X=fftn(([accX/mean(abs(accX)),accY/mean(abs(accY)),accZ/mean(abs(accZ))]));
    %%what we had^^
    %X = sum(fftn(([accX/mean(abs(accX)),accY/mean(abs(accY)),accZ/mean(abs(accZ))])),2);
    %X=fftn(sqrt(sum([ (accX/mean(abs(accX))).^2, (accY/mean(abs(accY))).^2 , (accZ/mean(abs(accZ))).^2 ],2)));
    %X=fftn(sum([accX/mean(abs(accX)),accY/mean(abs(accY)),accZ/mean(abs(accZ))],2));
    X=fftn(sum([accX,accY,accZ,mag(:,1),mag(:,2),mag(:,3)],2));
    %X=fft(sqrt(sum([accX.^2,accY.^2,accZ.^2],2)));
    %X=fft(accX.*accY.*accZ);
    n = length(accX);       % number of samples
    f = (0:n-1)*(fs/n);     % frequency range
    power = abs(X).^2/n;    % power of the DFT
    
    new=floor(length(f)/2);
    f=f(1:new);
    power=power(1:new);
    % ignore DC
    power(f<(.1))=0;
    %ignore over physical range
    power(f>(3))=0;
    
    f_where=f(power>=(1/4)*max(power));
    p_where=power(power>=(1/4)*max(power));
    f_dom_weighted=sum(f_where.*p_where/sum(p_where));
    f_dom=f(power==max(power));
    figure()
    plot(f, power);
    %pause()
    
    fprintf("file: %s \nf_dom: %d\n\n", data_files(i).name, f_dom);
    %close all
    fdoms=[fdoms,f_dom];

end



function [cps, int] = cycle_calculations(x,y,z,t)
    smallest_dt=.25;
    fs=1/(t(2)-t(1));
    %calculate cycles per second in t period
    [vals,locs] = findpeaks(x);
    xCycles=sum( vals>0 & [true; diff(locs)/fs>smallest_dt]);
    xlocs=locs( vals>0 & [true; diff(locs)/fs>smallest_dt]);
    [vals,locs] = findpeaks(y);
    yCycles=sum( vals>0 & [true; diff(locs)/fs>smallest_dt]);
    ylocs=locs( vals>0 & [true; diff(locs)/fs>smallest_dt]);
    [vals,locs] = findpeaks(z);
    zCycles=sum(vals>0 & [true; diff(locs)/fs>smallest_dt]);
    zlocs=locs(vals>0 & [true; diff(locs)/fs>smallest_dt]);
%     xCycles=sum(  [true; diff(locs)/fs>smallest_dt]);
%     xlocs=locs(  [true; diff(locs)/fs>smallest_dt]);
%     [vals,locs] = findpeaks(y);
%     yCycles=sum(  [true; diff(locs)/fs>smallest_dt]);
%     ylocs=locs(  [true; diff(locs)/fs>smallest_dt]);
%     [vals,locs] = findpeaks(z);
%     zCycles=sum( [true; diff(locs)/fs>smallest_dt]);
%     zlocs=locs( [true; diff(locs)/fs>smallest_dt]);
    
    cps=mean([xCycles,yCycles,zCycles])/t(end);
    %calculate pushing area in last n seconds
    %aproximate: numerical integration over radial frame. find center point
    %for each cycle and integrate change in r by change in theta
    
    %don't do first
    int=0;
    for i=2:length(xlocs)-1
        xhere=x(xlocs(i):xlocs(i+1));
        yhere=y(xlocs(i):xlocs(i+1));
        zhere=z(xlocs(i):xlocs(i+1));
        here=[xhere,yhere,zhere];
        %close all
        %plot3(here(:,1),here(:,2),here(:,3))
        
        if false
            center=mean(here);
            radii=here-center;
            dists=sqrt(   (radii(:,1).^2)   +    (radii(:,2).^2)  +    (radii(:,3).^2));
            speed=diff(radii,1);
            swept = cumsum(sqrt(speed(:,1).^2+speed(:,2).^2+speed(:,3).^2));
            angles = ((length(swept)-1)/length(swept))*2*pi*(swept)'/max(swept); %((length(swept)-1)/length(swept))
            %make it a closed figure
            dists(length(dists)+1)=dists(1);
            angles=[0,angles,2*pi];
            int(i)=calc_polar_area(dists, angles);
        end
        %% projection method
        if true
            %https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
            %fit plane to data ax+by+c=z    
            A=[xhere,yhere,ones(length(yhere),1)];
            B=zhere;
            coeffs = inv(A'*A)*A'*B;
            plane_origin = [0,0,coeffs(3)]; %c because x=y=0, []
            %here
            % project onto 2-space system
            %https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
            unit_normal = coeffs/norm(coeffs);
    %         hold on
    %         plot3([plane_origin(1),plane_origin(1)+unit_normal(1)],[plane_origin(2),plane_origin(2)+unit_normal(2)],[plane_origin(3),plane_origin(3)+unit_normal(3)],'--')
    %         hold off
            v = here-plane_origin;
            dist = v*unit_normal;
            projected_points = here - dist*unit_normal';
            here=projected_points;

            center=mean(here);

    %         figure()
    %         plot3(here(:,1),here(:,2),here(:,3))
    %         hold on
    %         scatter3(center(1),center(2),center(3))

            radii=here-center;
            dists=sqrt(   (radii(:,1).^2)   +    (radii(:,2).^2)  +    (radii(:,3).^2))';
            dists(length(dists)+1)=dists(1);
            speed=diff(radii,1);
            swept = cumsum(sqrt(speed(:,1).^2+speed(:,2).^2+speed(:,3).^2));
            angles = ((length(swept)-1)/length(swept))*2*pi*(swept)'/max(swept); %((length(swept)-1)/length(swept))
            angles=[0,angles,2*pi];
            %angles=2*pi*(1:length(dists))/length(dists); %this is an approximation- figure this out. something to do with speed vector
            %angles=angles(1:length(dists));%
            int(i)=calc_polar_area(dists, angles);
        end
    end
end
  
function [x_plane,y_plane,z_plane] = flatten_cloud(xhere,yhere,zhere)
    here=[xhere,yhere,zhere];
    A=[xhere,yhere,ones(length(yhere),1)];
    B=zhere;
    coeffs = inv(A'*A)*A'*B;
    plane_origin = [0,0,coeffs(3)]; %c because x=y=0, []
    %here
    % project onto 2-space system
    %https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
    unit_normal = coeffs/norm(coeffs);
    v = here-plane_origin;
    dist = v*unit_normal;
    projected_points = here - dist*unit_normal';
    x_plane=projected_points(:,1);
    y_plane=projected_points(:,2);
    z_plane=projected_points(:,3);
    

end

function int = calc_polar_area(dists, angles)
    %dists(length(dists)+1)=dists(1);
    int=0;
    %plot(dists.*cos(angles),dists.*sin(angles))
    for i=1:length(angles)-1
        dtheta=angles(i+1)-angles(i);
        h=(mean([dists(i),dists(i+1)]));
        dtheta(dtheta>2*pi)=dtheta(dtheta>2*pi)-2*pi;
        int=int+ .5*h*h*sin(dtheta);  
        %int=int+ .5*h*h*(dtheta);            
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
