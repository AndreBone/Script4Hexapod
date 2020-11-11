clear all
close all

%% settings
%range of values for r0, r1
r0_val =  1975;
r1_val =  1700;
%value of the angle between legs
sp0_val =3;
sp1_val =3.5:.1:4;
%range of values for the total height of the hexapod
h_val=1000:1:3000;

%Center of Mass in cyllindrical coords
zcm_val=-787.35-500;
rcm =   44.6128+200;
spcm=  214.0819;

%Uncertainties for CM position
ezcm = 500;
ercm = 200;

%Platform Position
x_val=0;
y_val=0;
z_val=0;
%Platform Rotation
th_val = 0;
ph_val = 0;
ps_val = 0;

%Platform positioning precision
dx=0.1; dy=0.1; dz=0.1; dth = 1/60; dph = 1/60; dps=1/60;

%Earthquake requirements
F_val = 3.6*[[1;0.3;0.3] [1;0.3;-0.3] [1;-0.3;0.3] [1;-0.3;-0.3] [-1;0.3;0.3] [-1;0.3;-0.3] [-1;-0.3;0.3] [-1;-0.3;-0.3]...
    [0.3;1;0.3] [0.3;1;-0.3] [-0.3;1;0.3] [-0.3;1;-0.3] [0.3;-1;0.3] [0.3;-1;-0.3] [-0.3;-1;0.3] [-0.3;-1;-0.3]...
    [0.3;0.3;1] [0.3;-0.3;1] [-0.3;0.3;1] [-0.3;-0.3;1] [0.3;0.3;-1] [0.3;-0.3;-1] [-0.3;0.3;-1] [-0.3;-0.3;-1]] + [0;0;-1];

%number of cycles to be run
count=size(r0_val,2)*size(r1_val,2)*size(sp0_val,2)*size(sp1_val,2)*size(h_val,2);
["Number of cycles to be run will be" num2str(count)]

config = nan(count,9);
MaxT_stg = zeros(size(zcm_val,2),6);

%% cycle

%geometry counter
j=0;

for r0 = r0_val
    for sp0 = sp0_val
        %Base Fixation Points
        p01 = r0*[cosd(30- sp0); sind(30- sp0); 0];
        p02 = r0*[cosd(30+ sp0); sind(30+ sp0); 0];
        p03 = r0*[cosd(150-sp0); sind(150-sp0); 0];
        p04 = r0*[cosd(150+sp0); sind(150+sp0); 0];
        p05 = r0*[cosd(270-sp0); sind(270-sp0); 0];
        p06 = r0*[cosd(270+sp0); sind(270+sp0); 0];
        for r1 = r1_val
            for sp1  =sp1_val
                for h = h_val
                    
                    if sp1==60-sp0
                        %This condition leads to singular matrices
                        continue
                    end
                    
                    %geometry counter
                    j=j+1;
                    
                    %platform Fixation Points
                    p11 = r1*[cosd(330+sp1); sind(330+sp1); 0];
                    p12 = r1*[cosd(90- sp1); sind(90- sp1); 0];
                    p13 = r1*[cosd(90+ sp1); sind(90+ sp1); 0];
                    p14 = r1*[cosd(210-sp1); sind(210-sp1); 0];
                    p15 = r1*[cosd(210+sp1); sind(210+sp1); 0];
                    p16 = r1*[cosd(330-sp1); sind(330-sp1); 0];
                    
                    for x = x_val
                        for y = y_val
                            for z = z_val+h
                                for th = th_val
                                    for ph = ph_val
                                        for ps = ps_val
                                            
                                            %Rotation matrices
                                            Rx = [1 0 0; 0 cosd(th) -sind(th); 0 sind(th) cosd(th)];
                                            Ry = [cosd(ph) 0 sind(ph); 0 1 0; -sind(ph) 0 cosd(ph)];
                                            Rz = [cosd(ps) -sind(ps) 0; sind(ps) cosd(ps) 0; 0 0 1];
                                            
                                            %Transformation matrix
                                            Tr = [Rx*Ry*Rz [x; y; z]; 0 0 0 1];
                                            
                                            %Applying Transformation matrix
                                            np11 = Tr*[p11; 1]; np11 = np11(1:3);
                                            np12 = Tr*[p12; 1]; np12 = np12(1:3);
                                            np13 = Tr*[p13; 1]; np13 = np13(1:3);
                                            np14 = Tr*[p14; 1]; np14 = np14(1:3);
                                            np15 = Tr*[p15; 1]; np15 = np15(1:3);
                                            np16 = Tr*[p16; 1]; np16 = np16(1:3);
                                            
                                            %Arm length
                                            v01 = np11-p01;
                                            v02 = np12-p02;
                                            v03 = np13-p03;
                                            v04 = np14-p04;
                                            v05 = np15-p05;
                                            v06 = np16-p06;
                                            
                                            %                     %Draw figure
                                            %                     figure(1)
                                            %                     fill3([p01(1); p02(1); p03(1); p04(1); p05(1); p06(1)], [p01(2); p02(2); p03(2); p04(2); p05(2); p06(2)], [p01(3); p02(3); p03(3); p04(3); p05(3); p06(3)],[ 0 ,148,255]/255)
                                            %                     hold on
                                            %                     fill3([p11(1); p12(1); p13(1); p14(1); p15(1); p16(1)], [p11(2); p12(2); p13(2); p14(2); p15(2); p16(2)], [p11(3); p12(3); p13(3); p14(3); p15(3); p16(3)],[123, 17, 66]/255)
                                            %
                                            %                     line([p01(1),p11(1)],[p01(2),p11(2)],[p01(3),p11(3)])
                                            %                     line([p02(1),p12(1)],[p02(2),p12(2)],[p02(3),p12(3)])
                                            %                     line([p03(1),p13(1)],[p03(2),p13(2)],[p03(3),p13(3)])
                                            %                     line([p04(1),p14(1)],[p04(2),p14(2)],[p04(3),p14(3)])
                                            %                     line([p05(1),p15(1)],[p05(2),p15(2)],[p05(3),p15(3)])
                                            %                     line([p06(1),p16(1)],[p06(2),p16(2)],[p06(3),p16(3)])
                                            %                     hold off
                                            %                     axis([-1850,1850,-1850,1850,0,2375])
                                            %                     daspect([1 1 1])
                                            %                     drawnow
                                            %                     print(['RH-Loop-' num2str(j)],'-dpng')
                                            
                                            
                                            
                                            %Force Analysis
                                            
                                            T_min = 10000;
                                            for F=F_val
                                                i=0;
                                                for zcm =zcm_val
                                                    i=i+1;
                                                    
                                                    %Center of mass position
                                                    pcm = Tr*[rcm*cosd(spcm); rcm*sind(spcm); zcm; 1]; pcm = pcm(1:3);
                                                    
                                                    %Now we'll solve the system MatrixA(VA) * Vector of Tensions in x (VTx) =
                                                    %Vector Solution (VS)
                                                    
                                                    VS=[-F;-cross(F,[x;y;z]-pcm)];
                                                    MA=[v01 v02 v03 v04 v05 v06;
                                                        cross(v01,[x;y;z]-np11) cross(v02,[x;y;z]-np12) cross(v03,[x;y;z]-np13) ...
                                                        cross(v04,[x;y;z]-np14) cross(v05,[x;y;z]-np15) cross(v06,[x;y;z]-np16)];
                                                    Tmat_stg    = [v01 v02 v03 v04 v05 v06]*diag(MA\VS);
                                                    
                                                    %Save the maximum force
                                                    %made by all of the
                                                    %arms, for all the
                                                    %forces in F_val
                                                    MaxT_stg(i,:) = max(MaxT_stg(i,:),[norm(Tmat_stg(:,1)), norm(Tmat_stg(:,2)), norm(Tmat_stg(:,3)),...
                                                        norm(Tmat_stg(:,4)), norm(Tmat_stg(:,5)), norm(Tmat_stg(:,6))]);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    %Filters the maximum force made by the arms for a
                    %certain zcm
                    MaxT = max(MaxT_stg,[],2);
                    
                    %chooses the zcm that minimizes the force
                    if min(MaxT) < T_min
                        T_min = min(MaxT);
                        spmin = spcm;
                        zcm_min = mean(zcm_val(MaxT==T_min));
                    end
                    
                    %Cler MaxT_stg
                    MaxT_stg = zeros(size(zcm_val,2),6);
                    
                    %Precision analysis
                    temp_der = 1e5;
                    
                    %platform Fixation Points
                    p11 = [r1*cosd(330+sp1); r1*sind(330+sp1); h];
                    p12 = [r1*cosd(90- sp1); r1*sind(90- sp1); h];
                    p13 = [r1*cosd(90+ sp1); r1*sind(90+ sp1); h];
                    p14 = [r1*cosd(210-sp1); r1*sind(210-sp1); h];
                    p15 = [r1*cosd(210+sp1); r1*sind(210+sp1); h];
                    p16 = [r1*cosd(330-sp1); r1*sind(330-sp1); h];
                    
                    %Arm length
                    v01 = p11-p01;
                    v02 = p12-p02;
                    v03 = p13-p03;
                    v04 = p14-p04;
                    v05 = p15-p05;
                    v06 = p16-p06;
                    
                    for drv=[[0; 0; dz; 0; 0; 0] [dx; 0; 0; 0; 0; 0] [0; dy; 0; 0; 0; 0] [0; 0; 0; dth; 0; 0] [0; 0; 0; 0; dph; 0] [0; 0; 0; 0; 0; dps]...
                            [-dx; 0; 0; 0; 0; 0] [0; -dy; 0; 0; 0; 0] [0; 0; -dz; 0; 0; 0] [0; 0; 0; -dth; 0; 0] [0; 0; 0; 0; -dph; 0] [0; 0; 0; 0; 0; -dps]]
                        Rx = [1 0 0; 0 cosd(drv(4)) -sind(drv(4)); 0 sind(drv(4)) cosd(drv(4))];
                        Ry = [cosd(drv(5)) 0 sind(drv(5)); 0 1 0; -sind(drv(5)) 0 cosd(drv(5))];
                        Rz = [cosd(drv(6)) -sind(drv(6)) 0; sind(drv(6)) cosd(drv(6)) 0; 0 0 1];
                        
                        %Transformation Matrix
                        Tr = [Rx*Ry*Rz [drv(1);drv(2);drv(3)]; 0 0 0 1];
                        
                        np11 = Tr*[p11; 1]; np11 = np11(1:3);
                        np12 = Tr*[p12; 1]; np12 = np12(1:3);
                        np13 = Tr*[p13; 1]; np13 = np13(1:3);
                        np14 = Tr*[p14; 1]; np14 = np14(1:3);
                        np15 = Tr*[p15; 1]; np15 = np15(1:3);
                        np16 = Tr*[p16; 1]; np16 = np16(1:3);
                        
                        nv01 = np11-p01;
                        nv02 = np12-p02;
                        nv03 = np13-p03;
                        nv04 = np14-p04;
                        nv05 = np15-p05;
                        nv06 = np16-p06;
                        
                        %Compute the difference
                        Jac = ([norm(nv01) norm(nv02) norm(nv03) norm(nv04) norm(nv05) norm(nv06)]-[norm(v01) norm(v02) norm(v03) norm(v04) norm(v05) norm(v06)]);
                        %Filter some small values
                        Jac(abs(Jac)<1e-10) = 0;
                        temp_der = min([temp_der abs(Jac(abs(Jac)>1e-10))]);
                        %temp_der = min([temp_der norm(Jac)]);
                    end
                    config(j,:) = [r0/1000, r1/1000, sp0, sp1, h/1000, zcm_min/1000, T_min, temp_der*1000, T_min];
                end
            end
            waitbar(j/count)
        end
    end
end

%% recalc
['r0 (m)' ' r1 (m)' ' sp0 (degrees)' ' sp1 (degrees)' ' h (m)' ' zcm_min (m) ' , ' load ', ' precision (um)', ' Cost Function']
uconfig=config(config(:,end)==min(config(:,end)),:)

r0  =uconfig(1)*1000;
r1  =uconfig(2)*1000;
sp0 =uconfig(3);
sp1 =uconfig(4);
h  =uconfig(5)*1000;

p01 = r0*[cosd(30- sp0); sind(30- sp0); 0];
p02 = r0*[cosd(30+ sp0); sind(30+ sp0); 0];
p03 = r0*[cosd(150-sp0); sind(150-sp0); 0];
p04 = r0*[cosd(150+sp0); sind(150+sp0); 0];
p05 = r0*[cosd(270-sp0); sind(270-sp0); 0];
p06 = r0*[cosd(270+sp0); sind(270+sp0); 0];

p11 = [r1*cosd(330+sp1); r1*sind(330+sp1); h];
p12 = [r1*cosd(90- sp1); r1*sind(90- sp1); h];
p13 = [r1*cosd(90+ sp1); r1*sind(90+ sp1); h];
p14 = [r1*cosd(210-sp1); r1*sind(210-sp1); h];
p15 = [r1*cosd(210+sp1); r1*sind(210+sp1); h];
p16 = [r1*cosd(330-sp1); r1*sind(330-sp1); h];

v01 = p11-p01;
v02 = p12-p02;
v03 = p13-p03;
v04 = p14-p04;
v05 = p15-p05;
v06 = p16-p06;

figure(1)
fill3([p01(1); p02(1); p03(1); p04(1); p05(1); p06(1)], [p01(2); p02(2); p03(2); p04(2); p05(2); p06(2)], [p01(3); p02(3); p03(3); p04(3); p05(3); p06(3)],[ 0 ,148,255]/255)
hold on
fill3([p11(1); p12(1); p13(1); p14(1); p15(1); p16(1)], [p11(2); p12(2); p13(2); p14(2); p15(2); p16(2)], [p11(3); p12(3); p13(3); p14(3); p15(3); p16(3)],[123, 17, 66]/255)

line([p01(1),p11(1)],[p01(2),p11(2)],[p01(3),p11(3)])
line([p02(1),p12(1)],[p02(2),p12(2)],[p02(3),p12(3)])
line([p03(1),p13(1)],[p03(2),p13(2)],[p03(3),p13(3)])
line([p04(1),p14(1)],[p04(2),p14(2)],[p04(3),p14(3)])
line([p05(1),p15(1)],[p05(2),p15(2)],[p05(3),p15(3)])
line([p06(1),p16(1)],[p06(2),p16(2)],[p06(3),p16(3)])
hold off
%axis([-1850,1850,-1850,1850,0,2375])
daspect([1 1 1])
drawnow

clear('MaxT')
MaxT_stg = zeros(201,6);
MaxT_Com = 0;
MaxT_Str = 0;
j=0;
FV_C=0;
FV_S=0;
for rcmt = [max(rcm-ercm,0):5:rcm+ercm]
    j=j+1;
    for F=F_val
        i=0;
        for zcmt =[zcm-ezcm:5:zcm+ezcm]
            i=i+1;
            
            pcm = [rcmt*cosd(spcm); rcmt*sind(spcm); zcmt+h];
            
            VS=[-F;-cross(F,[x;y;z]-pcm)];
            MA=[v01 v02 v03 v04 v05 v06;
                cross(v01,[x;y;z]-p11) cross(v02,[x;y;z]-p12) cross(v03,[x;y;z]-p13) ...
                cross(v04,[x;y;z]-p14) cross(v05,[x;y;z]-p15) cross(v06,[x;y;z]-p16)];
            Tmat_stg    = [v01 v02 v03 v04 v05 v06]*diag(MA\VS);
            MaxT_stg(i,:) = max(MaxT_stg(i,:),[norm(Tmat_stg(:,1)), norm(Tmat_stg(:,2)), norm(Tmat_stg(:,3)),...
                norm(Tmat_stg(:,4)), norm(Tmat_stg(:,5)), norm(Tmat_stg(:,6))]);
            MaxT_Com = max([MaxT_Com, norm(Tmat_stg(:,1))*(dot(v01,Tmat_stg(:,1))>0), norm(Tmat_stg(:,2))*(dot(v02,Tmat_stg(:,2))>0),...
                norm(Tmat_stg(:,3))*(dot(v03,Tmat_stg(:,3))>0), norm(Tmat_stg(:,4))*(dot(v04,Tmat_stg(:,4))>0),...
                norm(Tmat_stg(:,5))*(dot(v05,Tmat_stg(:,5))>0), norm(Tmat_stg(:,6))*(dot(v06,Tmat_stg(:,6))>0)]);
            MaxT_Str = max([MaxT_Str, norm(Tmat_stg(:,1))*(dot(v01,Tmat_stg(:,1))<0), norm(Tmat_stg(:,2))*(dot(v02,Tmat_stg(:,2))<0),...
                norm(Tmat_stg(:,3))*(dot(v03,Tmat_stg(:,3))<0), norm(Tmat_stg(:,4))*(dot(v04,Tmat_stg(:,4))<0),...
                norm(Tmat_stg(:,5))*(dot(v05,Tmat_stg(:,5))<0), norm(Tmat_stg(:,6))*(dot(v06,Tmat_stg(:,6))<0)]);
        end
    end
    MaxT(j,:) = max(MaxT_stg,[],2)';
end

figure(2)
imagesc(zcm-ezcm:5:zcm+ezcm,max(rcm-ercm,0):5:rcm+ercm,MaxT)
drawnow
xlabel('z deviation (mm)')
ylabel('Lateral deviation (mm)')
cb = colorbar;
title(cb,'Load (METIS weight)')
title('Arm load under earthquake')

["Under earthquake" min(min(MaxT)) max(max(MaxT))]

clear('MaxT')
MaxT_stg = zeros(201,6);
j=0;
for rcmt = [max(rcm-ercm,0):5:rcm+ercm]
    j=j+1;
    F=[0;0;-1];
    i=0;
    for zcmt =[zcm-ezcm:5:zcm+ezcm]
        i=i+1;
        
        pcm = [rcmt*cosd(spcm); rcmt*sind(spcm); zcmt+h];
        
        VS=[-F;-cross(F,[x;y;z]-pcm)];
        MA=[v01 v02 v03 v04 v05 v06;
            cross(v01,[x;y;z]-p11) cross(v02,[x;y;z]-p12) cross(v03,[x;y;z]-p13) ...
            cross(v04,[x;y;z]-p14) cross(v05,[x;y;z]-p15) cross(v06,[x;y;z]-p16)];
        Tmat_stg    = [v01 v02 v03 v04 v05 v06]*diag(MA\VS);
        MaxT_stg(i,:) = max(MaxT_stg(i,:),[norm(Tmat_stg(:,1)), norm(Tmat_stg(:,2)), norm(Tmat_stg(:,3)),...
            norm(Tmat_stg(:,4)), norm(Tmat_stg(:,5)), norm(Tmat_stg(:,6))]);
    end
    MaxT(j,:) = max(MaxT_stg,[],2)';
end

figure(3)
imagesc(zcm-ezcm:5:zcm+ezcm,max(rcm-ercm,0):5:rcm+ercm,MaxT)
drawnow
xlabel('z deviation (mm)')
ylabel('Lateral deviation (mm)')
cb = colorbar;
title(cb,'Load (METIS weight)')
title('Arm load under normal gravity')

["Under normal gravity" min(min(MaxT)) max(max(MaxT))]

%% positions

lmin=min([norm(v01), norm(v02), norm(v03), norm(v04), norm(v05), norm(v06)]);
lmax=max([norm(v01), norm(v02), norm(v03), norm(v04), norm(v05), norm(v06)]);
vamax=0;
hamax=0;


for pos = [[-50;0;0;0;0;0], [50;0;0;0;0;0], [0;-50;0;0;0;0], [0;50;0;0;0;0], [0;0;-50;0;0;0], [0;0;400;0;0;0],... 
           [0;0;0;2;0;0], [0;0;0;-2;0;0], [0;0;0;0;2;0], [0;0;0;0;-2;0], [0;0;0;0;0;2], [0;0;0;0;0;-2]]
    x=pos(1);
    y=pos(2);
    z=pos(3);
    th=pos(4);
    ph=pos(5);
    ps=pos(6);
    
    %Transformation Matrix
    Rx = [1 0 0; 0 cosd(th) -sind(th); 0 sind(th) cosd(th)];
    Ry = [cosd(ph) 0 sind(ph); 0 1 0; -sind(ph) 0 cosd(ph)];
    Rz = [cosd(ps) -sind(ps) 0; sind(ps) cosd(ps) 0; 0 0 1];
    
    Tr = [Rx*Ry*Rz [x;y;z]; 0 0 0 1];
    
    np11 = Tr*[p11; 1]; np11 = np11(1:3);
    np12 = Tr*[p12; 1]; np12 = np12(1:3);
    np13 = Tr*[p13; 1]; np13 = np13(1:3);
    np14 = Tr*[p14; 1]; np14 = np14(1:3);
    np15 = Tr*[p15; 1]; np15 = np15(1:3);
    np16 = Tr*[p16; 1]; np16 = np16(1:3);
    
    nv01 = np11-p01;
    nv02 = np12-p02;
    nv03 = np13-p03;
    nv04 = np14-p04;
    nv05 = np15-p05;
    nv06 = np16-p06;
    
    lmin=min([lmin norm(nv01), norm(nv02), norm(nv03), norm(nv04), norm(nv05), norm(nv06)]);
    lmax=max([lmax norm(nv01), norm(nv02), norm(nv03), norm(nv04), norm(nv05), norm(nv06)]);
    
    for v=[[v01;nv01] [v02;nv02] [v03;nv03] [v04;nv04] [v05;nv05] [v06;nv06]]
        a=v(1:3)/norm(v(1:3));
        b=[a(1); a(2); -(a(1)^2+a(2)^2)/a(3)];
        b=b/norm(b);
        c=[-a(2)/sqrt(a(1)^2+a(2)^2); a(1)/sqrt(a(1)^2+a(2)^2); 0];
        
        d=[c b a]\v(4:6);
        d=d/norm(d);
        %d=[sin(ha); -sin(va)cos(ha); cos(va) cos(ha)]
        vamax=max([vamax, acosd(sqrt(d(2)^2+d(3)^2))]);
        hamax=max([hamax, abs(atand(d(2)/d(3)))]);
    end
end


["length range" lmin lmax lmax-lmin]
["angles of joints" vamax hamax]
