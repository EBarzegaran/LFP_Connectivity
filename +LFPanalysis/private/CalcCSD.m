function csddata = CalcBipolar(data)
% Calculates bipolar data for the input matrix

% INPUT: 
    % data: nchan x nsample x ntrial matrix

% Author: Elham Barzegaran, September, 2019
%% calcualte the bipolar data
data1 = zeros(size(data,1)+2,size(data,2),size(data,3));
data2 = data1;
data3 = data1;

data1(3:end,:,:) = data;% shift down
data1(2,:,:) = data(1,:,:);

data2(1:end-2,:,:) = data;% shift up
data2(end-1,:,:) = data(end,:,:);

data3(2:end-1,:,:)=data;% original


csddata = data2+data1-2*data3;
csddata = csddata(2:end-1,:,:);
end