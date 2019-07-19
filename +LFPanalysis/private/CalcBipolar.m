function bdata = CalcBipolar(data)
% Calculates bipolar data for the input matrix

% INPUT: 
    % data: nchan x nsample x ntrial matrix

% Author: Elham Barzegaran, July, 2019
%% calcualte the bipolar data
data1 = zeros(size(data,1)+2,size(data,2),size(data,3));
data2 = data1;

data1(3:end,:,:) = data;% shift down
data1(2,:,:) = data(1,:,:);

data2(1:end-2,:,:) = data;% shift up
data2(end-1,:,:) = data(end,:,:);

bdata = data2-data1;
bdata = bdata(2:end-1,:,:);
end