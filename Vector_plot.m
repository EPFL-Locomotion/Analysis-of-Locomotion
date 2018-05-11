%
Vector=zeros(4*size(S6_FLOAT.T_01.Raw.Kin.LHIP,1),3);

k=1;

for i=1:size(S6_FLOAT.T_01.Raw.Kin.LHIP)
Vector(k,:)=S6_FLOAT.T_01.Raw.Kin.LHIP(i,:);
Vector(k+1,:)=S6_FLOAT.T_01.Raw.Kin.LKNE(i,:);
Vector(k+2,:)=S6_FLOAT.T_01.Raw.Kin.LANK(i,:);
Vector(k+3,:)=S6_FLOAT.T_01.Raw.Kin.LTOE(i,:);

k=k+4;
end

k=1;
for i=1:4:size(S6_FLOAT.T_01.Raw.Kin.LHIP,1)/10
    figure;
    
    plot3(Vector(i:3+i,1),Vector(i:i+3,2),Vector(i:i+3,3));
    hold on
    scatter3(S6_FLOAT.T_01.Raw.Kin.LHIP(k,1),S6_FLOAT.T_01.Raw.Kin.LHIP(k,2),S6_FLOAT.T_01.Raw.Kin.LHIP(k,3));
    hold on
    scatter3(S6_FLOAT.T_01.Raw.Kin.LKNE(k,1),S6_FLOAT.T_01.Raw.Kin.LKNE(k,2),S6_FLOAT.T_01.Raw.Kin.LKNE(k,3))
    hold on
    scatter3(S6_FLOAT.T_01.Raw.Kin.LANK(k,1),S6_FLOAT.T_01.Raw.Kin.LANK(k,2),S6_FLOAT.T_01.Raw.Kin.LANK(k,3))
    hold on
    scatter3(S6_FLOAT.T_01.Raw.Kin.LTOE(k,1),S6_FLOAT.T_01.Raw.Kin.LTOE(k,2),S6_FLOAT.T_01.Raw.Kin.LTOE(k,3))
    k=k+1;
end