clear; clc; close all;

[f v n]=stlread('Bronze_statue_STL.stl');
a=f.Points;
X=a(:,1);
Y=a(:,2);
Z=a(:,3);

patch('faces',f.ConnectivityList,'vertices',f.Points,'edgecolor','none','facecolor','red');
% patch('faces',f.ConnectivityList,'vertices',f.Points,'edgecolor','white','facecolor','red');
camlight

% return



% 평면 정하기 %

the1=rand*pi/2;                                     % the1 각
phi1=rand*pi/2;                                     % phi1 각

b1=[cos(the1)*cos(phi1)
    sin(the1)*cos(phi1)
    sin(phi1)];                                     % b1 벡터      
bx1=[-sin(the1)*cos(phi1)-cos(the1)*sin(phi1)
    cos(the1)*cos(phi1)-sin(the1)*sin(phi1)
    cos(phi1)];                                     % bx1 벡터

by1=cross(b1,bx1);                                  % bx2 벡터

the2=rand*pi/2;                                     % the2 각
phi2=rand*pi/2;                                     % phi2 각

b2=[cos(the2)*cos(phi2)
    sin(the2)*cos(phi2)
    sin(phi2)];                                     % b2 벡터
bx2=[-sin(the2)*cos(phi2)-cos(the2)*sin(phi2)
    cos(the2)*cos(phi2)-sin(the2)*sin(phi2)
    cos(phi2)];                                     % bx1 벡터

by2=cross(b2,bx2);                                  % bx2 벡터

real=[X Y Z];
figure
subplot(2,2,[1,2]);

mDataX1=real*bx1/sum(bx1.^2);                   % bx1 벡터에 투사한 좌표
mDataY1=real*by1/sum(by1.^2);                   % by1 벡터에 투사한 좌표
% plot(mDataX1,mDataY1,'.b');                      % bx1 벡터와 by1벡터의 좌표로 그리기

mConnect1=delaunay(mDataX1,mDataY1);
triplot(mConnect1,mDataX1,mDataY1);

hold on
title(sprintf('%1.0f %1.0f',the1*180/pi,phi1*180/pi));


subplot(2,2,[3,4])
% figure

mDataX2=real*bx2/sum(bx2.^2);                   % bx2 벡터에 투사한 좌표      
mDataY2=real*by2/sum(by2.^2);                   % by2 벡터에 투사한 좌표
% plot(mDataX2,mDataY2,'.b');                      % bx2 벡터와 by2 벡터의 좌표로 그리기
mConnect2=delaunay(mDataX2,mDataY2);
triplot(mConnect2,mDataX2,mDataY2);
hold on
title(sprintf('%1.0f %1.0f',the2*180/pi,phi2*180/pi));

% return


figure

B=[bx1 by1 bx2 by2];
b=[sum(bx1.^2)*mDataX1  sum(by1.^2)*mDataY1 ...
    sum(bx2.^2)*mDataX2  sum(by2.^2)*mDataY2];
make=B'\b';
% plot3(make(1,:),make(2,:),make(3,:),'.r')        % 투사한 데이터 그리기
% hold on

% box on
% grid on
make=make';
patch('faces',f.ConnectivityList,'vertices',make,'edgecolor','none','facecolor','blue');
% % patch('faces',f.ConnectivityList,'vertices',f.Points,'edgecolor','white','facecolor','red');
camlight
