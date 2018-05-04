% Segment the data into gait cycles 

LTOEGait = NO_FLOAT_CRUTCHES.T_02.Raw.Kin.LTOE(:,3);

% Scan from beginning to foot off (from data) for each foot off
% Initial contact (IC) = Foot Strike

[LeftIC(1,2),LeftIC(1,1)]  = max(LTOEGait(1:300));
[LeftIC(2,2),LeftIC(2,1)]  = max(LTOEGait(301:582));
[LeftIC(3,2),LeftIC(3,1)]  = max(LTOEGait(583:883));
[LeftIC(4,2),LeftIC(4,1)]  = max(LTOEGait(884:1182));
[LeftIC(5,2),LeftIC(5,1)]  = max(LTOEGait(1183:1472));
[LeftIC(6,2),LeftIC(6,1)]  = max(LTOEGait(1473:1775));
[LeftIC(7,2),LeftIC(7,1)]  = max(LTOEGait(1776:2085));
[LeftIC(8,2),LeftIC(8,1)]  = max(LTOEGait(2086:2155));


%correction for time mismatch (beginning of each scan)

LeftIC(2,1)  = LeftIC(2,1) + 300;
LeftIC(3,1)  = LeftIC(3,1) + 582;
LeftIC(4,1)  = LeftIC(4,1) + 883;
LeftIC(5,1)  = LeftIC(5,1) + 1182;
LeftIC(6,1)  = LeftIC(6,1) + 1472;
LeftIC(7,1)  = LeftIC(7,1) + 1775;
LeftIC(8,1)  = LeftIC(8,1) + 2085;

% Terminal contact (IC) = Foot Off

[LeftTC(1,2),LeftTC(1,1)]  = min(LTOEGait(1:300));
[LeftTC(2,2),LeftTC(2,1)]  = min(LTOEGait(301:590));
[LeftTC(3,2),LeftTC(3,1)]  = min(LTOEGait(591:900));
[LeftTC(4,2),LeftTC(4,1)]  = min(LTOEGait(901:1200));
[LeftTC(5,2),LeftTC(5,1)]  = min(LTOEGait(1201:1500));
[LeftTC(6,2),LeftTC(6,1)]  = min(LTOEGait(1501:1800));
[LeftTC(7,2),LeftTC(7,1)]  = min(LTOEGait(1801:2200));

%correction for time mismatch (beginning of each scan)

LeftTC(2,1)  = LeftTC(2,1) + 300;
LeftTC(3,1)  = LeftTC(3,1) + 591;
LeftTC(4,1)  = LeftTC(4,1) + 900;
LeftTC(5,1)  = LeftTC(5,1) + 1200;
LeftTC(6,1)  = LeftTC(6,1) + 1500;
LeftTC(7,1)  = LeftTC(7,1) + 1800;

hold on

plot(NO_FLOAT_CRUTCHES.T_02.Raw.Kin.LTOE(:,3));
title('T02 Left Toe Angular momentum')
legend('z axis');
xlabel('Time(s)');
ylabel('Angular velocity');
plot(LeftIC(1:8, 1), LeftIC(1:8, 2),'or');
plot(LeftTC(1:7, 1), LeftTC(1:7, 2), 'og');

hold off

disp(LeftIC);
disp(LeftTC);

