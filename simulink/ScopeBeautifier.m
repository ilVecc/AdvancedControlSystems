%% prettify pose data

var = out.PoseData;

n = 3;
if size(var.signals, 2) > 3
    m = 2;
else
    m = 1;
end

idx = reshape(1:(n*m), m, n)';
for i=1:length(var.signals)
    subplot(n,m,idx(i));
    plot(var.time, squeeze(var.signals(i).values(1,1,:)),'b');
    hold on;
    plot(var.time, squeeze(var.signals(i).values(1,2,:)),'r');
    title(var.signals(i).label)
end
sgtitle('Pose') 

clear n m i idx var


%% prettify force data

var = out.ForceData;

n = 3;
if size(var.signals, 2) > 3
    m = 2;
else
    m = 1;
end

idx = reshape(1:(n*m), m, n)';
for i=1:length(var.signals)
    subplot(n,m,idx(i));
    plot(var.time, squeeze(var.signals(i).values(:,1)),'m');
    title(var.signals(i).label)
end
sgtitle('Force') 

clear n m i idx var


%% prettify force (no torque) data

var = out.ForceData;

n = 3;
if size(var.signals, 2) > 3
    m = 2;
else
    m = 1;
end

idx = reshape(1:(n*m), m, n)';
for i=1:length(var.signals)
    subplot(n,m,idx(i));
    plot(var.time, squeeze(var.signals(i).values(1,1,:)),'m');
    hold on;
    plot(var.time, squeeze(var.signals(i).values(2,1,:)),'g');
    title(var.signals(i).label)
end

clear n m i idx var
