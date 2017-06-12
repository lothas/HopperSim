classdef HopperModel
    % Version 0.2 - 14/09/2012
    
    % 6 DOF monoped hopper
    % 2 DOF position + Orientation of body and leg
    % + length of leg + length of spring
    
    % 14/09
    % Applied a torque to the foot when the terrain is inclined
    
    properties(Constant)
        % Link Identifiers
        L_body=1;
        L_leg1=2;
        L_leg2=3;
        
        % Extremities Identifiers
        E_leg=4;
        E_foot=5;
        
        COM=6;
        
        % Gravity
        g=9.81;
    end
    
    properties
        mBody=5;        % Body link mass
        lBody=0.6;      % Body link length
        IBody=0;        % Body link inertia
        mLeg1=2;        % Leg upper link mass
        lLeg1=0.25;     % Leg upper link length
        ILeg1=0;        % Leg upper link inertia
        mLeg2=2;        % Leg lower link mass
        lLeg2=0.25;     % Leg lower link length
        ILeg2=0;        % Leg lower link inertia
        Spr_l0=0.2;     % Spring initial length
        Spr_k=1000;     % Spring constant
        Damping=5;      % Spring damper
                
        % Ground interaction
        Floor=[];       % Terrain object
        OnGround=0;     % set to 1 if foot is on the ground
        PosFoot0=[0;0];  % position of foot on the ground
        Foot_Radius=0.1; % Radius of the foot point
                
        % Inputs
        Torque=0;
        Force=0;
                
        % Render parameters
        RenderObj;
        
        LinkWidth=0.04;
        Colors={[0.9,0.3,0.2],[0.7,0.7,0.8],[0.3, 0.3, 0.4],[0.5, 0.5, 1]};
                
        CircRes=6;
        LinkRes=5;
        
        LineWidth=1;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [HM] = HopperModel(varargin)
            switch nargin
                case 0
                    HM=HM.SetInertia();
                case 1
                case 2
                case 3
                case 4
            end     
        end
        
        function [HM] = SetInertia(HM)
            HM.IBody=HM.mBody*HM.lBody^2/12;
            HM.ILeg1=HM.mLeg1*HM.lLeg1^2/12;
            HM.ILeg2=HM.mLeg2*HM.lLeg2^2/12;
        end
        
        function [xy] = GetPos(HM,q,which)
            switch which
                % %%%%%%%%%% Links %%%%%%%%%% %
                case HM.L_body
                    xy=[q(1)+(q(4)+q(5))*sin(q(3));
                        q(2)+(q(4)+q(5))*cos(q(3))];
                case HM.L_leg1
                    xy=[q(1)+(q(4)+HM.lLeg1/2)*sin(q(3));
                        q(2)+(q(4)+HM.lLeg1/2)*cos(q(3))];
                case HM.L_leg2
                    xy=[q(1)+(q(4)+q(5)-HM.lLeg2/2)*sin(q(3));
                        q(2)+(q(4)+q(5)-HM.lLeg2/2)*cos(q(3))];

                % %%%%%%% Extremities %%%%%% %
                case HM.E_leg
                    xy=[q(1)+q(4)*sin(q(3));
                        q(2)+q(4)*cos(q(3))];
                case HM.E_foot
                    if HM.OnGround
                        xy=HM.PosFoot0;
                    else
                        xy=[q(1); q(2)];
                    end
                    
                % %%%%%% Center of Mass %%%%%% %
                case HM.COM 
                    [xy1] = HM.GetPos(q,HM.L_body);
                    [xy2] = HM.GetPos(q,HM.L_leg1);
                    [xy3] = HM.GetPos(q,HM.L_leg2);
                    xy=[(xy1(1)*HM.mBody+xy2(1)*HM.mLeg1+xy3(1)*HM.mLeg2)/(HM.mBody+HM.mLeg1+HM.mLeg2);
                        (xy1(2)*HM.mBody+xy2(2)*HM.mLeg1+xy3(2)*HM.mLeg2)/(HM.mBody+HM.mLeg1+HM.mLeg2)];
            end
        end
        
        function [xy] = GetVel(HM,q,which)
            switch which
                % %%%%%%%%%% Links %%%%%%%%%% %
                case HM.L_body
                    xy=[q(6+1)+(q(6+4)+q(6+5))*sin(q(3))+(q(4)+q(5))*cos(q(3))*q(6+3);
                        q(6+2)+(q(6+4)+q(6+5))*cos(q(3))-(q(4)+q(5))*sin(q(3))*q(6+3)];
                case HM.L_leg1
                    xy=[q(6+1)+q(6+4)*sin(q(3))+(q(4)+HM.lLeg1/2)*cos(q(3))*q(6+3);
                        q(6+2)+q(6+4)*cos(q(3))-(q(4)+HM.lLeg1/2)*sin(q(3))*q(6+3)];
                case HM.L_leg2
                    xy=[q(6+1)+(q(6+4)+q(6+5))*sin(q(3))+(q(4)+q(5)-HM.lLeg2/2)*cos(q(3))*q(6+3);
                        q(6+2)+(q(6+4)+q(6+5))*cos(q(3))-(q(4)+q(5)-HM.lLeg2/2)*sin(q(3))*q(6+3)];

                % %%%%%%% Extremities %%%%%% %
                case HM.E_leg
                    xy=[q(6+1)+q(6+4)*sin(q(3))+q(4)*cos(q(3))*q(6+3);
                        q(6+2)+q(6+4)*cos(q(3))-q(4)*sin(q(3))*q(6+3)];
                case HM.E_foot
                    if HM.OnGround
                        xy=[0; 0];
                    else
                        xy=[q(6+1); q(6+2)];
                    end
                    
                % %%%%%% Center of Mass %%%%%% %
                case HM.COM 
                    xy=[0; 0];
            end
        end
                
        function [qdot]=Derivative4D(HM,t,q) %#ok<INUSL>
            % q = [ theta1, lS, lT, theta2 ] 
            [H,h,Q]=HM.GetMatrices4D(q);
            
            q_dot2=pinv(H)*(Q-h);
            q_dot2(3)=HM.Force;
            
            qdot=[q(5:8); q_dot2];
        end
        
        function [qdot]=Derivative5D(HM,t,q) %#ok<INUSL>
            % q = [ x, y, theta1, lT, theta2 ] 
            [H,h,Q]=HM.GetMatrices5D(q);
            
            q_dot2=pinv(H)*(Q-h);
            q_dot2(4)=HM.Force;
            
            qdot=[q(6:10); q_dot2];
        end
        
        function [H,h,Q]=GetMatrices4D(HM,q)
            % %%%%%%%%%%%%%%%%%%%%% MATRIX H %%%%%%%%%%%%%%%%%%%%% %
            H=zeros(4,4);
            % [row] is for the equation
            % [col] is for the coordinate [theta1, lS, lT, theta2]
            
            
            % theta1 equation
            H(1,1)=HM.IBody+HM.ILeg1+HM.ILeg2+q(3)^2*(HM.mBody+HM.mLeg2)+(HM.lLeg1^2*HM.mLeg1)/4+(HM.lLeg2^2*HM.mLeg2)/4+...
                q(2)*(HM.lLeg1*HM.mLeg1-HM.lLeg2*HM.mLeg2)+q(2)^2*(HM.mBody+HM.mLeg1+HM.mLeg2)+...
                q(2)*q(3)*(2*HM.mBody+2*HM.mLeg2)-HM.lLeg2*HM.mLeg2*q(3);
            H(1,4)=HM.IBody;

            % lS equation
            H(2,2)=HM.mBody+HM.mLeg1+HM.mLeg2;
            H(2,3)=HM.mBody+HM.mLeg2;
            
            % lT equation
            H(3,2)=H(2,3);
            H(3,3)=H(2,3);
            
            % theta2 equation
            H(4,1)=H(1,4);
            H(4,4)=HM.IBody;
                        
            % %%%%%%%%%%%%%%%%%%%%% VECTOR h %%%%%%%%%%%%%%%%%%%%% %
            h=zeros(4,1);
            
            kk1=HM.mBody+HM.mLeg1+HM.mLeg2;
            kk2=HM.mBody+HM.mLeg2;
            kk3=HM.lLeg1*HM.mLeg1/2-HM.lLeg2*HM.mLeg2/2;
            
            h(1)=q(4+1)*(q(2)*(q(4+2)*2*kk1+q(4+3)*2*kk2)+q(4+2)*2*kk3+q(3)*(q(4+2)*2*kk2+q(4+3)*2*kk2)-HM.lLeg2*HM.mLeg2*q(4+3))-...
                HM.g*sin(q(1))*(q(3)*kk2+kk3+q(2)*kk1)+HM.Spr_k*(HM.Spr_l0-q(2))*HM.Foot_Radius*sin(HM.Floor.SurfSlope(HM.PosFoot0(1))+q(1));
            h(2)=(-kk3-kk1*q(2)-kk2*q(3))*q(4+1)^2-HM.Spr_k*(HM.Spr_l0-q(2))+HM.g*cos(q(1))*kk1+HM.Damping*q(4+2);
            h(3)=-q(4+1)^2*(-(HM.lLeg2*HM.mLeg2)/2+kk2*(q(3)+q(2)))+HM.g*cos(q(1))*kk2;
            h(4)=0;
            
            % %%%%%%%%%%%%%%%%%%%%% VECTOR F %%%%%%%%%%%%%%%%%%%%% %
            Q=zeros(4,1);
            
            % theta2 equation
            Q(4)=HM.Torque;
            
            % lT equation
            Q(3)=HM.Force;
        end
        
        function [H,h,Q]=GetMatrices5D(HM,q)
            % %%%%%%%%%%%%%%%%%%%%% MATRIX H %%%%%%%%%%%%%%%%%%%%% %
            H=zeros(5,5);
            % [row] is for the equation
            % [col] is for the coordinate [x, y, theta1, lT, theta2]
            
            kk1=HM.mBody+HM.mLeg1+HM.mLeg2;
            kk2=HM.mBody+HM.mLeg2;
            kk3=HM.lLeg1*HM.mLeg1/2-HM.lLeg2*HM.mLeg2/2;
            
            % x equation
            H(1,1)=kk1;
            H(1,3)=cos(q(3))*(HM.Spr_l0*kk1+kk3+kk2*q(4));
            H(1,4)=sin(q(3))*kk2;
            
            % y equation
            H(2,2)=kk1;
            H(2,3)=-sin(q(3))*(HM.Spr_l0*kk1+kk3+kk2*q(4));
            H(2,4)=cos(q(3))*kk2;
            
            % theta1 equation
            H(3,1)=H(1,3);
            H(3,2)=H(2,3);
            H(3,3)=HM.IBody+HM.ILeg1+HM.ILeg2+HM.Spr_l0^2*kk1+(HM.lLeg1^2*HM.mLeg1)/4+(HM.lLeg2^2*HM.mLeg2)/4+kk2*q(4)^2+HM.Spr_l0*2*kk3+2*HM.Spr_l0*q(4)*kk2-HM.lLeg2*HM.mLeg2*q(4);
            H(3,5)=HM.IBody;
                        
            % lT equation
            H(4,1)=H(1,4);
            H(4,2)=H(2,4);
            H(4,4)=kk2;
            
            % theta2 equation
            H(5,3)=H(3,5);
            H(5,5)=HM.IBody;
            
            % %%%%%%%%%%%%%%%%%%%%% VECTOR h %%%%%%%%%%%%%%%%%%%%% %
            h=zeros(5,1);
            
            h(1)=2*q(5+3)*q(5+4)*cos(q(3))*kk2-q(5+3)^2*sin(q(3))*(HM.Spr_l0*kk1+kk3+kk2*q(4));
            h(2)=-cos(q(3))*(HM.Spr_l0*kk1+kk3+kk2*q(4))*q(5+3)^2-2*q(5+4)*sin(q(3))*kk2*q(5+3)+HM.g*kk1;
            h(3)=q(5+3)*q(5+4)*(2*HM.Spr_l0*kk2-HM.lLeg2*HM.mLeg2+2*kk2*q(4))-HM.g*sin(q(3))*(HM.Spr_l0*kk1+kk3+kk2*q(4));
            h(4)=HM.g*cos(q(3))*kk2-q(5+3)^2*(HM.Spr_l0*kk2-(HM.lLeg2*HM.mLeg2)/2+kk2*q(4));
            h(5)=0;
            
            % %%%%%%%%%%%%%%%%%%%%% VECTOR F %%%%%%%%%%%%%%%%%%%%% %
            Q=zeros(5,1);
            
            % lT equation
            Q(4)=HM.Force;
            
            % theta2 equation
            Q(5)=HM.Torque;
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value isterminal direction] = Events(HM, q)
            value=ones(2,1);
            isterminal=ones(2,1);
            direction=ones(2,1);
            
            PosFoot=HM.GetPos(q,HM.E_foot);
        
            % Check for foot contact
            value(1)=PosFoot(2)-HM.Floor.Surf(PosFoot(1));
            
            % Check for foot detachment
            value(2)=HM.Spr_l0-HM.GetSpringLength(q);
        end
        
        function [HM] = Render(HM,q)
            % Get all the positions required
            % Body
            PosBody=HM.GetPos(q,HM.L_body);
            % Leg upper link
            PosLegL1=HM.GetPos(q,HM.L_leg1);
            % Leg lower link
            PosLegL2=HM.GetPos(q,HM.L_leg2);
            % Spring start
            PosSpringS=HM.GetPos(q,HM.E_leg);
            % Spring end
            PosSpringE=HM.GetPos(q,HM.E_foot);
            
            if isempty(HM.RenderObj)
                HM.RenderObj.Body=HM.DrawLink(PosBody, q(3)+q(6), 3, HM.lBody, 1.1*HM.LinkWidth, HM.Colors{1});
                HM.RenderObj.Leg1=HM.DrawLink(PosLegL1, q(3)+pi/2, 2, HM.lLeg1, 0.7*HM.LinkWidth, HM.Colors{2});
                HM.RenderObj.Leg2=HM.DrawLink(PosLegL2, q(3)+pi/2, 1, HM.lLeg2, HM.LinkWidth, HM.Colors{3});
                hold on
                HM.RenderObj.Spring=HM.DrawSpring(PosSpringS, PosSpringE, []);
                axis equal
            else
                if ishandle(HM.RenderObj.Body.Geom)==0
                    % Graphic handles were destroyed and need to be renewed
                    HM.RenderObj=[];
                    HM=HM.Render(q);
                else
                    HM.DrawLink(PosBody, q(3)+q(6), HM.RenderObj.Body);
                    HM.DrawLink(PosLegL1, q(3)+pi/2, HM.RenderObj.Leg1);
                    HM.DrawLink(PosLegL2, q(3)+pi/2, HM.RenderObj.Leg2);
                    HM.DrawSpring(PosSpringS, PosSpringE, HM.RenderObj.Spring);
                end
            end            
        end
       
        function [ res ] = DrawLink(HM, Pos, Or, z_or_Obj, Length, Width, Color)
            switch nargin
                case 4
                    Txy=makehgtform('translate',[Pos(1) Pos(2) 0]);
                    Rz=makehgtform('zrotate',-Or);
                    set(z_or_Obj.Trans,'Matrix',Txy*Rz);
                    res=1;
                case 7
                    res.Trans=hgtransform('Parent',gca);
                    Txy=makehgtform('translate',[Pos(1) Pos(2) z_or_Obj]);
                    Rz=makehgtform('zrotate',-Or);

                    coordX=[-Length/2, -Length/2, Length/2, Length/2];
                    coordY=[-Width/2, Width/2, Width/2, -Width/2];
                    coordZ=[z_or_Obj, z_or_Obj, z_or_Obj, z_or_Obj];

                    res.Geom=patch(coordX,coordY,coordZ,Color);
                    set(res.Geom,'EdgeColor',[0 0 0]);
                    set(res.Geom,'LineWidth',2*HM.LineWidth);

                    set(res.Geom,'Parent',res.Trans);
                    set(res.Trans,'Matrix',Txy*Rz);
            end
        end
        
        function [ res ] = DrawSpring(HM, PosS, PosE, Obj)
            NumTurns=5;
            
            Center=(PosS+PosE)/2;
            Length=sqrt((PosE(1)-PosS(1))^2+(PosE(2)-PosS(2))^2);
            Orientation=atan2(PosE(2)-PosS(2),PosE(1)-PosS(1));

            coordX=zeros(1,NumTurns*2+2);
            coordY=zeros(1,NumTurns*2+2);

            Step=Length/NumTurns;

            coordX(1)=-Length/2;
            coordX(end)=Length/2;
            for i=1:NumTurns
                coordX(2*i)=coordX(1)+Step*(1/4+i-1);
                coordX(2*i+1)=coordX(1)+Step*(3/4+i-1);
                coordY(2*i)=HM.LinkWidth/2;
                coordY(2*i+1)=-HM.LinkWidth/2;
            end
            
            % Rotate spring
            Coords=[coordX; coordY];
            CoordsFoot=Coords(:,end);
            CoordsFoot(1)=CoordsFoot(1)-HM.LinkWidth/2;
            RotMatrix=[cos(Orientation) -sin(Orientation);
                       sin(Orientation)  cos(Orientation)];
            Coords=repmat(Center,1,NumTurns*2+2)+RotMatrix*Coords;
            CoordsFoot=Center+RotMatrix*CoordsFoot;
            
            if isempty(Obj)
                res=zeros(2,1);
                res(1)=plot(Coords(1,:), Coords(2,:), 'Color', [0, 0, 0], 'LineWidth', 2*HM.LineWidth);
                res(2)=plot(CoordsFoot(1), CoordsFoot(2), '.', 'Color', [0.3, 0.3, 0.35], 'MarkerSize', 600*HM.Foot_Radius);
            else
                set(Obj(1),'XData',Coords(1,:));
                set(Obj(1),'YData',Coords(2,:));
                set(Obj(2),'XData',CoordsFoot(1));
                set(Obj(2),'YData',CoordsFoot(2));
                res=1;
            end
        end
    end
    
    methods(Static)
        function [y] = Ground(x)
            y=-0.02+0.05*x;
        end
        
        function [SR] = SatRamp(x)
            SR=min(max(100*x,0),1);
        end
    end
end