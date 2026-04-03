%% 读取TPXO10数据，并写入boundary文件中
clc
clear

%% 读取TPXO潮位数据
filename = 'H:\BaiduSyncdisk\05Delft\MIKE_Python\Astronomical\TPXO10_atlas_v2_nc\grid_tpxo10atlas_v2.nc';
folder = 'H:\BaiduSyncdisk\05Delft\MIKE_Python\Astronomical\TPXO10_atlas_v2_nc';

S = dir(fullfile(folder, '*.nc'));   % 结构体数组
names = {S.name};                    % 仅文件名（cell）

bounds=[105 125;15 25]';
%[Z,x_or_lon,y_or_lat] = tmd_data(filename,variable)

x_or_lon = double(ncread(filename,'lon_z'));
y_or_lat = double(ncread(filename,'lat_z'));

dx = abs(diff(x_or_lon(1:2)));
dy = abs(diff(y_or_lat(1:2)));

% Get xlimits (xl) and ylimits (yl) of input coordinates + buffer:
xl = [min(bounds(:,1))-2*dx max(bounds(:,1))+2*dx];
yl = [min(bounds(:,2))-2*dy max(bounds(:,2))+2*dy];

% Region of rows and columns of pixels to read:
ri=find((y_or_lat>=yl(1))&(y_or_lat<=yl(2)));
ci=find((x_or_lon>=xl(1))&(x_or_lon<=xl(2)));

x_or_lon = x_or_lon(ci);
y_or_lat = y_or_lat(ri);

conList= {'m2','s2','n2','k2','k1','o1','p1','q1','mf','mm','m4','ms4','mn4','2n2','s1'};
NCons = length(conList);
Z = NaN(numel(ri),numel(ci),NCons);
P = NaN(numel(ri),numel(ci),NCons);

for k = 1:NCons
    placInd = k;
    filename= fullfile(folder, names{contains(names,['_',conList{k}])+contains(names,'h')==2});%根据variable检索文件序号
    Z(:,:,placInd) = abs(complex(double(permute(ncread(filename,'hRe',[ri(1) ci(1)],[numel(ri) numel(ci)]),[1 2])),...
        double(permute(ncread(filename,'hIm',[ri(1) ci(1)],[numel(ri) numel(ci)]),[1,2]))));%ncread(source,varname,start,count)
    P(:,:,placInd) = angle(complex(double(permute(ncread(filename,'hRe',[ri(1) ci(1)],[numel(ri) numel(ci)]),[1 2])),...
        double(permute(ncread(filename,'hIm',[ri(1) ci(1)],[numel(ri) numel(ci)]),[1 2]))));
end
clearvars -except P Z x_or_lon y_or_lat conList NCons
%% Interpolate
% 读取需要生成的边界条件点
folderPli = 'H:\BaiduSyncdisk\Phd2\MIKE_Python\Project11_1d2d.dsproj_data\FlowFM\input';
plif = dir(fullfile(folderPli, '*.pli'));   % 结构体数组
namesPli = {plif.name};                    % 仅文件名（cell）
fid = fopen('WaterLevel.bc','a');   % 覆盖写
for i=1:length(namesPli) %分别给每个pli文件写入bc
    data = readmatrix(fullfile(folderPli,namesPli{i}),FileType="text");
    data=data(:,1:2);
    xii=data(:,1);
    yii=data(:,2);
    % 加密点
    ds=11000;% unit=m
    Pnew = data(1,:);%第一个点不变
    for m = 1:size(data,1)-1
        Pm = data(m,:); Pn = data(m+1,:);
        v  = Pn - Pm;
        L  = hypot(v(1), v(2));
        if L < ds %点的间距小于指定间距
            sNew=L;
        else %如果间距大于指点间距则加密
            sNew = (0:ds:L)'; %新的线段端点
            sNew(end+1,1)=L; %保证包含最后一个点,把第一个点去掉，避免重复
            sNew(1,:)=[];
        end
        t=sNew/L; %每个线段占几分之几
        Pnew =[Pnew;Pm + t .* v]; % 新的点 
    end
    [xi,yi]=fromXian2WGS(Pnew(:,1),Pnew(:,2));%读取x和y并转换成经纬度,便于插值
    
    zi = nan(size(xi,1),NCons);% amplitude
    phi = nan(size(xi,1),NCons); %phase
    for k = 1:NCons
        zi(:,k) = interp2(x_or_lon,y_or_lat,Z(:,:,k),xi,yi);
        phi(:,k) = interp2(x_or_lon,y_or_lat,P(:,:,k),xi,yi);
    end
    zi=zi*0.001;%从mm转到m
    phi=phi*180/pi;%弧度转到度
    % 把经纬度坐标转回Xian
    [xp,yp]=fromWGS2Xian(xi,yi);

    %% 写入pli文件
    fid1 = fopen(['Boundary',sprintf('%02d', i),'.pli'],'a');   % 覆盖写
    LineStart=1;
    np=length(Pnew);%该线段的点数
    fprintf(fid1, '%s\n', ['Boundary',sprintf('%02d', i)]);
    fprintf(fid1, '%s\t%s\n', num2str(np),    '2');
    % 追加写入矩阵
    for p=1:np
        fprintf(fid1, '%s\t%s\t%s\n',xp(p), yp(p),['Boundary',sprintf('%02d', i),'_',sprintf('%04d', p)]);
    end
    fclose(fid1);
    %% 写入boundary文件
    for np=1:length(xi)
        fprintf(fid, '%s\n','[forcing]');
        fprintf(fid, '%s\n', ['Name         =Boundary',sprintf('%02d', i),'_',sprintf('%04d', np)]);
        fprintf(fid, '%s\n', 'Function = astronomic');
        fprintf(fid, '%s\n', 'Quantity = astronomic component');
        fprintf(fid, '%s\n', 'Unit = -');
        fprintf(fid, '%s\n', 'Quantity = waterlevelbnd amplitude');
        fprintf(fid, '%s\n', 'Unit = m');
        fprintf(fid, '%s\n', 'Quantity = waterlevelbnd phase');
        fprintf(fid, '%s\n', 'Unit = deg');
        for ncon=1:NCons
            fprintf(fid, '%s\t%s\t%s\n',upper(conList{ncon}), num2str(zi(np,ncon)),num2str(phi(np,ncon)));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
