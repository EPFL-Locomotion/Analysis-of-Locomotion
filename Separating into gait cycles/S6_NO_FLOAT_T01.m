% Segment the data into gait cycles 

LTOEGait = S6_NO_FLOAT.T_01.Raw.Kin.LTOE(:,3);

% Scan from beginning to foot off (from data) for each foot off
% Initial contact (IC) = Foot Strike

[LeftIC(1,2),LeftIC(1,1)]  = max(LTOEGait(1:200));
[LeftIC(2,2),LeftIC(2,1)]  = max(LTOEGait(201:320));
[LeftIC(3,2),LeftIC(3,1)]  = max(LTOEGait(321:450));
[LeftIC(4,2),LeftIC(4,1)]  = max(LTOEGait(451:560));


%correction for time mismatch (beginning of each scan)

LeftIC(2,1)  = LeftIC(2,1) + 200;
LeftIC(3,1)  = LeftIC(3,1) + 320;
LeftIC(4,1)  = LeftIC(4,1) + 450;

% Terminal contact (IC) = Foot Off

[LeftTC(1,2),LeftTC(1,1)]  = min(LTOEGait(1:300));
[LeftTC(2,2),LeftTC(2,1)]  = min(LTOEGait(301:450));
[LeftTC(3,2),LeftTC(3,1)]  = min(LTOEGait(451:550));
[LeftTC(4,2),LeftTC(4,1)]  = min(LTOEGait(550:650));


%correction for time mismatch (beginning of each scan)

LeftTC(2,1)  = LeftTC(2,1) + 300;
LeftTC(3,1)  = LeftTC(3,1) + 450;
LeftTC(4,1)  = LeftTC(4,1) + 550;

hold on

plot(S6_NO_FLOAT.T_01.Raw.Kin.LTOE(:,3));
title('T01 Left Toe Angular momentum')
legend('z axis');
xlabel('Time(s)');
ylabel('Angular velocity');
plot(LeftIC(1:4, 1), LeftIC(1:4, 2),'or');
plot(LeftTC(1:4, 1), LeftTC(1:4, 2), 'og');

hold off

disp(LeftIC(1:4));
disp(LeftTC(1:4));

