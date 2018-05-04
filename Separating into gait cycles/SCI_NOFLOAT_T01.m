% Segment the data into gait cycles 

LTOEGait = NO_FLOAT_CRUTCHES.T_01.Raw.Kin.LTOE(:,3);

% Scan from beginning to foot off (from data) for each foot off
% Initial contact (IC) = Foot Strike

[LeftIC(1,2),LeftIC(1,1)]  = max(LTOEGait(1:233));
[LeftIC(2,2),LeftIC(2,1)]  = max(LTOEGait(234:508));
[LeftIC(3,2),LeftIC(3,1)]  = max(LTOEGait(509:792));
[LeftIC(4,2),LeftIC(4,1)]  = max(LTOEGait(793:1009));
[LeftIC(5,2),LeftIC(5,1)]  = max(LTOEGait(1010:1390));
[LeftIC(6,2),LeftIC(6,1)]  = max(LTOEGait(1390:1688));
[LeftIC(7,2),LeftIC(7,1)]  = max(LTOEGait(1689:1963));
[LeftIC(8,2),LeftIC(8,1)]  = max(LTOEGait(1964:2284));
[LeftIC(9,2),LeftIC(9,1)]  = max(LTOEGait(2285:2340));

%correction for time mismatch (beginning of each scan)

LeftIC(2,1)  = LeftIC(2,1) + 233;
LeftIC(3,1)  = LeftIC(3,1) + 508;
LeftIC(4,1)  = LeftIC(4,1) + 792;
LeftIC(5,1)  = LeftIC(5,1) + 1009;
LeftIC(6,1)  = LeftIC(6,1) + 1390;
LeftIC(7,1)  = LeftIC(7,1) + 1688;
LeftIC(8,1)  = LeftIC(8,1) + 1963;
LeftIC(9,1)  = LeftIC(9,1) + 2284;

% Terminal contact (IC) = Foot Off

[LeftTC(1,2),LeftTC(1,1)]  = min(LTOEGait(1:233));
[LeftTC(2,2),LeftTC(2,1)]  = min(LTOEGait(234:508));
[LeftTC(3,2),LeftTC(3,1)]  = min(LTOEGait(509:792));
[LeftTC(4,2),LeftTC(4,1)]  = min(LTOEGait(793:1050));
[LeftTC(5,2),LeftTC(5,1)]  = min(LTOEGait(1051:1390));
[LeftTC(6,2),LeftTC(6,1)]  = min(LTOEGait(1390:1688));
[LeftTC(7,2),LeftTC(7,1)]  = min(LTOEGait(1689:1963));
[LeftTC(8,2),LeftTC(8,1)]  = min(LTOEGait(1964:2284));

%correction for time mismatch (beginning of each scan)

LeftTC(2,1)  = LeftTC(2,1) + 233;
LeftTC(3,1)  = LeftTC(3,1) + 508;
LeftTC(4,1)  = LeftTC(4,1) + 792;
LeftTC(5,1)  = LeftTC(5,1) + 1050;
LeftTC(6,1)  = LeftTC(6,1) + 1390;
LeftTC(7,1)  = LeftTC(7,1) + 1688;
LeftTC(8,1)  = LeftTC(8,1) + 1963;

hold on

plot(NO_FLOAT_CRUTCHES.T_01.Raw.Kin.LTOE(:,3));
title('T01 Left Toe Angular momentum')
legend('z axis');
xlabel('Time(s)');
ylabel('Angular velocity');
plot(LeftIC(1:end, 1), LeftIC(1:end, 2),'or');
plot(LeftTC(1:end, 1), LeftTC(1:end, 2), 'og');

hold off

disp(LeftIC);
disp(LeftTC);

