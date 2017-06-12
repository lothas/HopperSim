function [  ] = Simulation(  )
	% Version 0.1 - 11/09/2012
    
Video=0; % Output to video

%% ----------------- Parameter definitions -----------------
% Load robot model
Robot=HopperModel();

% Set initial conditions
InitCond=[0,0.3,0,0.2,0.3,0,...
          0.,0,0.,0,0,0];
Robot.OnGround=0;

Floor = Terrain(1,0.1,2);
% Floor.SetSmoothness(200);
Robot.Floor=Floor;

% Indexes within state vector X for each component
IdRobot=1:10;

% Set simulation time vector
tstart=0;
tend=80;
tstep=0.02;
tspan=tstart:tstep:tend;

% Set value for manual stop of simulation (1 = stop)
StopSim=0;

% Render parameters
Graphics=1;
Follow=1;

scrsz=0;
FigWin=0;
tCOM=0;
COM0=0;
TimeDisp=0;
Once=1;

FlMin0=-0.75;
FlMax0=0.75;
HeightMin0=-0.1;
HeightMax0=1.1;
FlMin=FlMin0;
FlMax=FlMax0;
HeightMin=HeightMin0;
HeightMax=HeightMax0;

%% ----------------- Functions definitions ----------------- 

    % Dynamical equations for the whole system
    function [Xdot] = Derivative(t, X)
        Xdot=zeros(12,1);
        [Robot.Force,Robot.Torque] = Controller(t, X);
        
        if Robot.OnGround
            X4D=[X(3:6);X(9:12)];
            X4Ddot=Robot.Derivative4D(t, X4D);
            Xdot(1:2)=0;
            Xdot(3:6)=X4Ddot(1:4);
            Xdot(7:8)=0;
            Xdot(9:12)=X4Ddot(5:8);
        else
            X5D=[X(1:3);X(5:9);X(11:12)];
            X5Ddot=Robot.Derivative5D(t, X5D);
            Xdot(1:3)=X5Ddot(1:3);
            Xdot(4)=0;
            Xdot(5:9)=X5Ddot(4:8);
            Xdot(10)=0;
            Xdot(11:12)=X5Ddot(9:10);
        end
    end

ForSw=Robot.OnGround;
PistPos0=0.3;
PistPos1=0.35;
PistP=500;
PistD=50;
LegP=250;
LegD=30;
LegAngle=0.01;
FwdSpeed=0.1;
tGo=0;
delta_tGo=0.03;
fbr=ones(5,1);

    function [Force,Torque] = Controller(t, X)
        % Add up to 10% noise to "measured" signal
        if t<20
            FwdSpeed=0.1*sin(2*pi*t/20);
        else
            FwdSpeed=0;
        end       
        
        if ForSw==1 || t<tGo
            Force=PistP*(PistPos1-X(5)*fbr(1))-PistD*X(6+5)*fbr(2);
            Torque=-LegP*(X(6)+X(3))*fbr(3)-LegD*2*(X(6+6)+X(3+6))*fbr(4);
        else
            [xy]=Robot.GetVel(X,Robot.L_body);
            LegAngle=xy(1)/2*fbr(1)-FwdSpeed;
            Force=PistP*(PistPos0-X(5)*fbr(2))-PistD*X(6+5)*fbr(3);
            Torque=LegP*(X(3)*fbr(4)+LegAngle)+LegD*X(6+3)*fbr(5);
        end
%         Torque=0;
    end
        

    % Event function for the whole system
    function [value isterminal direction] = EventFunc(t, X) %#ok<INUSL>
        numEvents=1+1; % +1 event to stop simulation on the spot
        value=ones(numEvents,1);
        isterminal=ones(numEvents,1);
        direction=-ones(numEvents,1);
        
        % Check for foot contact / detachment
        [xy]=Robot.GetPos(X,Robot.E_leg);
%         if xy(2)>Robot.Spr_l0 && t>0.4
%             d=1;
%         end
        if Robot.OnGround
            value(1)=Robot.Spr_l0-X(4);
        else
            value(1)=X(2)-Floor.Surf(xy(1));
        end
        
        if StopSim==1
            value(numEvents) = 0;
            direction(numEvents) = 0;
        end
    end

if Video==1
    aviobj = VideoWriter('SimVideo','Motion JPEG AVI');
    aviobj.FrameRate = 1/tstep/2; % 1/2 speed
    % Open file for writing
    open(aviobj);
end

    % Real time plot for the whole system
    function status = RealTimePlot(t,X,flag)
        if strcmp(flag,'done')
            % Finish simulation
            return
        end
               
        if Graphics
            if strcmp(flag,'init')
                if Once
                    % Initialize
                    scrsz = get(0, 'ScreenSize');
                    FigWin=figure();
                    set(FigWin,'Position', [150 100 scrsz(3)-500 scrsz(4)-300]);
                    set(gca,'LooseInset',get(gca,'TightInset')*2)

                    hold on
                    axis([FlMin FlMax HeightMin HeightMax]);
                    axis equal

                    [COM0]=Robot.GetPos(X,Robot.COM); 
                    COM0(2)=Floor.Surf(COM0(1));
                    tCOM = hgtransform('Parent',gca);

                    % Draw Robot
                    Robot=Robot.Render(X(IdRobot));
                    % Draw Floor
                    Floor=Floor.Render(FlMin,FlMax);

                    % Display time
%                     TimeDisp=text(0, -0.6*L, sprintf('t=%.2f',t),'HorizontalAlignment','center','Parent',tCOM);
                    TimeDisp=text(0.9*FlMax, 0.9, sprintf('t=%.2f',t),'HorizontalAlignment','center','Parent',tCOM);

                    uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                        'Position', [scrsz(3)-700 scrsz(4)-360 100 30],...
                        'Callback', @StopButtonCallback);

                    Once=0;
                end
                
                return
            end

            if Follow
                [COM]=Robot.GetPos(X,Robot.COM);
                COM(2)=max(Floor.Surf(COM(1)),COM(2)-1);
                FlMin=COM(1)+FlMin0;
                FlMax=COM(1)+FlMax0;
                HeightMin=COM(2)+HeightMin0;
                HeightMax=COM(2)+HeightMax0;
                
                TCOMx = makehgtform('translate',[COM(1)-COM0(1) COM(2)-COM0(2) 0]);
                set(tCOM,'Matrix',TCOMx);
                
                axis([FlMin FlMax HeightMin HeightMax]);
            end
                
            % Draw Robot
            Robot.Render(X(IdRobot));
            % Draw Floor
            Floor=Floor.Render(FlMin,FlMax);
   
            % Display time
%             if ischar(CurSpeed)
%                 set(TimeDisp,'String',sprintf(['t=%.2f - Osc.=%.3f\n',...
%                                                'Slope=%.2f\n',...
%                                                'Speed= %s'],t(1),X(5),Floor.SurfSlope(Robot.xS)*180/pi,CurSpeed));
%             else
%                 set(TimeDisp,'String',sprintf(['t=%.2f - Osc.=%.3f\n',...
%                                                'Slope=%.2f\n',...
%                                                'Speed= %.3f m/s (%.3f km/h)'],t(1),X(5),Floor.SurfSlope(Robot.xS)*180/pi,CurSpeed,CurSpeed*3.6));
%             end
            set(TimeDisp,'String',sprintf('t=%.2f',t(1)));
            
            drawnow
            
            if Video==1
                Frame=getframe;
                % Add frame to video
                writeVideo(aviobj,Frame);
            end
        end
        
        status=StopSim;
    end

    function StopButtonCallback(t,X) %#ok<INUSD>
        StopSim=1;
        close all
    end  % StopButtonCallback

%% ----------------- Simulation -----------------

options=odeset('MaxStep',tstep/10,'RelTol',.1e-7,'AbsTol',.1e-8, 'OutputFcn', @RealTimePlot, 'Events', @EventFunc);
[TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,InitCond,options); %#ok<ASGLU>

T=TTemp;
X=XTemp;
% Save torques
% TorquesTemp=zeros(2,length(TTemp));
% for j=1:length(TTemp)
%     TorquesTemp(:,j)=CPG.NeurOutput();
% end
% TorquesP=TorquesTemp;

while TTemp(end)<tspan(end-1) && StopSim==0
    for i=1:size(IE,1)
        switch IE(i)
            case 1 
                if Robot.OnGround % Foot detached
                    Robot.OnGround=0;
                    ForSw=0;
                    tGo=TTemp(end)+delta_tGo;
                    % Transform from 4D to 5D system
                    % Calculate velocity for foot after detachment
                    [xy]=Robot.GetVel(XTemp(end,:),Robot.L_body);
                    xtnew=xy(1)-(XTemp(end,4)+XTemp(end,5))*cos(XTemp(end,3))*XTemp(end,6+3);
                    ytnew=xy(2)+(XTemp(end,4)+XTemp(end,5))*sin(XTemp(end,3))*XTemp(end,6+3);
                    IC=[XTemp(end,1:6),xtnew,ytnew,XTemp(end,9:12)];
                else % Foot touched ground   
                    % Set new foothold
                    Robot.PosFoot0=Robot.GetPos(X(end,IdRobot),Robot.E_foot);
                    Robot.OnGround=1;
                    ForSw=1;
                    LegAngle=-LegAngle;
                    % Transform from 5D to 4D system
                    % velocity of foot will be transfered to spring
                    % compression speed
                    Speed=sqrt(XTemp(end,6+1)^2+XTemp(end,6+2)^2);
                    % and we also need to update the leg angular velocity
                    q3=XTemp(end,3);
                    q4=XTemp(end,4);
                    q5=XTemp(end,5);
                    q1tb=-XTemp(end,6+1);
                    q2tb=-XTemp(end,6+2);
                    q3tb=XTemp(end,6+3);
                    kk1=Robot.mBody+Robot.mLeg1+Robot.mLeg2;
                    kk2=Robot.mBody+Robot.mLeg2;
                    kk3=Robot.lLeg1*Robot.mLeg1-Robot.lLeg2*Robot.mLeg2;
                    kk4=(Robot.lLeg1^2*Robot.mLeg1)/4+(Robot.lLeg2^2*Robot.mLeg2)/4;
                    AngVel=(q3tb*(q5^2*kk2+kk4+q4*2*kk3+q4^2*kk1+q4*q5*2*kk2-Robot.lLeg2*Robot.mLeg2*q5)+...
                        q2tb*sin(q3)*(q5*kk2+kk3+q4*kk1)-q1tb*cos(q3)*(q5*kk2+kk3+q4*kk1))/...
                        (q5^2*kk2+kk4+q4*2*kk3+q4^2*kk1+q4*q5*2*kk2-Robot.lLeg2*Robot.mLeg2*q5);
                    IC=[XTemp(end,1:6),0,0,AngVel,-Speed,XTemp(end,11:12)];
                    % q = [ theta1, lS, lT, theta2 ] 
                    % q = [ x, y, theta1, lT, theta2 ] 
                end
            fbr=1+(0.5-rand(5,1))/10;
        end
    end
    
    
       
    % Check "End-Game" conditions
%     MaxAngle=pi/3+Floor.SurfSlope(Robot.xS);
%     if abs(Xnew(1))>MaxAngle
%         EndReached=1;
%         EndText='Stance leg angle too large';
%         break;
%     end
%     if abs(Xnew(2))>MaxAngle
%         EndReached=2;
%         EndText='Swing leg angle too large';
%         break;
%     end
%     if abs(Xnew(2)-Xnew(1))<0.00001
%         EndReached=3;
%         EndText='Step length too small';
%         break;
%     end
    
    % Continue simulation
    tspan=TTemp(end):tstep:tend;
    [TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,IC,options); %#ok<ASGLU>
    
    T=[T; TTemp]; %#ok<AGROW>
    X=[X; XTemp]; %#ok<AGROW>
    
    % Save torques
%     TorquesTemp=zeros(2,length(TTemp));
%     for j=1:length(TTemp)
%         TorquesTemp(:,j)=CPG.NeurOutput();
%     end
%     TorquesP=[TorquesP, TorquesTemp]; %#ok<AGROW>
end

% plot(T,X(:,19:21));

if Video==1
    % Close file
    close(aviobj);
end
end

