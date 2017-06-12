classdef Terrain
    % Version 0.21 - 15/07/2012
    
    % Different terrains: inclined plane, sinusoidal
    %                     infinite parabolla, finite parabolla
    
    properties
        Type=0;
        % 0 - inclined plane
        % 1 - sinusoidal
        % 2 - infinite parabolla
        % 3 - finite parabolla
        
        sinAmp=0.1;
        sinFreq=1;
        
        parK=0.025 % approx 1.5 deg/m
        
        incline=0;
        start_slope=0;
        end_slope=0;
        
        end_x=0;
        end_y=0;
        start_x=0;
        start_y=0;
        
        % Render parameters
        FloorStep=0.05;
        VertLines=10;
        FloorColor=[0.1,0.4,0];

        FloorLine=0;
        FloorVLx=zeros(1,10);
        FloorVL=zeros(1,10);
        
        LineWidth=1;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [Te] = Terrain(varargin)
            switch nargin
                case 0
                    Te; %#ok<VUNUS>
                case 1
                    Te.Type=varargin{1};
                case 2
                    Te.Type=varargin{1};
                    Te.start_slope=varargin{2};
                    Te.end_slope=varargin{2};
                case 3
                    Te.Type=varargin{1};
                    switch varargin{1}
                        case 1
                            Te.sinAmp=varargin{2};
                            Te.sinFreq=varargin{3};
                        case 2
                            Te.start_slope=varargin{2};
                            Te.parK=varargin{3};
                        case 3
                            Te.start_slope=varargin{2};
                            Te.end_slope=varargin{3};
                        otherwise
                            Te.start_slope=varargin{2};
                            Te.end_slope=varargin{2};
                    end
                case 4
                    Te.Type=varargin{1};
                    Te.start_slope=varargin{2};
                    Te.end_slope=varargin{3};
                    Te.parK=varargin{4};
            end     
            Te=Te.SetVertLines(Te.VertLines);
            Te=SetEndConditions(Te);
        end
        
        function [Te] = SetEndConditions(Te)
            if Te.Type==2
                Te.incline=sign(Te.parK);
                Te.parK=Te.parK/Te.incline;
            else
                Te.incline=sign(Te.end_slope-Te.start_slope);
            end
            
            if Te.incline==0
                Te.end_x=0;
                Te.end_y=0;
            else
                % end_x is actually the distance required
                % from x=0 to reach the desired angle
                Te.end_x=Te.start_x+(tan(Te.end_slope*pi/180)-tan(Te.start_slope*pi/180))/(Te.incline*Te.parK);
%                 Te.end_y=Te.start_y+(tan(Te.end_slope*pi/180)^2-tan(Te.start_slope*pi/180)^2)/2/Te.parK;
                Te.end_y=Te.start_y+tan(Te.start_slope*pi/180)*(Te.end_x-Te.start_x)+Te.incline*Te.parK/2*(Te.end_x-Te.start_x)^2;
            end
        end
        
        function [Te] = SetSmoothness(Te,Kslope)
            % The higher the value of Kslope, the faster the slope
            % will change from start_slope to end_slope
            % Default is 1.5 degrees per meter
            Te.parK=Kslope;
            Te=Te.SetEndConditions();
        end
        
        function [Te] = SetVertLines(Te,NLines)
            Te.VertLines=NLines;
            Te.FloorVLx=zeros(1,NLines);
            Te.FloorVL=zeros(1,NLines);
        end
        
        function [y, Trans] = Surf(Te,x)
            switch Te.Type
                case 0 % inclined plane
                    [y, Trans] = Te.Surf0(x);
                case 1 % sinusoidal
                    [y, Trans] = Te.Surf1(x);
                case 2 % infinite parabolla
                    [y, Trans] = Te.Surf2(x);
                case 3 % finite parabolla
                    [y, Trans] = Te.Surf3(x);
            end
        end

        function [alpha]=SurfSlope(Te,x)
            switch Te.Type
                case 0 % inclined plane
                    [alpha] = Te.SurfSlope0(x);
                case 1 % sinusoidal
                    [alpha] = Te.SurfSlope1(x);
                case 2 % infinite parabolla
                    [alpha] = Te.SurfSlope2(x);
                case 3 % finite parabolla
                    [alpha] = Te.SurfSlope3(x);
            end
        end
        
        % 0 - inclined plane
        function [y, Trans] = Surf0(Te,x)
            alpha = Te.SurfSlope(x);
            y=Te.start_y+x*tand(Te.end_slope);
            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
        end

        function [alpha]=SurfSlope0(Te,x) %#ok<INUSD>
            alpha=Te.end_slope*pi/180;
        end
        
        % 1 - sinusoidal
        function [y, Trans] = Surf1(Te,x)
            alpha = Te.SurfSlope(x);
            if x<Te.start_x
                y=0;            
            else
                y=Te.sinAmp*(1-cos(Te.sinFreq*(x-Te.start_x)));
            end

            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
        end

        function [alpha]=SurfSlope1(Te,x)
            if x<Te.start_x
                alpha=0;
            else
                alpha=Te.sinAmp*Te.sinFreq*sin(Te.sinFreq*(x));
            end
        end
        
        % 2 - infinite parabolla
        function [y, Trans] = Surf2(Te,x)
            alpha = Te.SurfSlope(x);
            if x<Te.start_x
                y=tan(alpha)*(x-Te.start_x);            
            else
                y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x)+Te.incline*Te.parK/2*(x-Te.start_x)^2;
            end

            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
        end

        function [alpha]=SurfSlope2(Te,x)
            if x<Te.start_x
                alpha=Te.start_slope*pi/180;
            else
                alpha=atan(tan(Te.start_slope*pi/180)+Te.incline*Te.parK*(x-Te.start_x));
            end
        end
        
        % 3 - finite parabolla
        function [y, Trans] = Surf3(Te,x)
            alpha = Te.SurfSlope(x);
            if x<Te.start_x
                y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x);            
            else
                if x<Te.end_x
                    y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x)+Te.incline*Te.parK/2*(x-Te.start_x)^2;
                else
                    y=Te.end_y+tan(Te.end_slope*pi/180)*(x-Te.end_x);
                end
            end

            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
        end

        function [alpha]=SurfSlope3(Te,x)
            if x<Te.start_x
                alpha=Te.start_slope*pi/180;
            else
                if x<Te.end_x
                    alpha=atan(tan(Te.start_slope*pi/180)+Te.incline*Te.parK*(x-Te.start_x));
                else
                    alpha=Te.end_slope*pi/180;
                end
            end
        end
        
        function [Te]=Render(Te,Min,Max)
            FloorX=Min:Te.FloorStep:Max;
            FN=length(FloorX);
            FloorY=zeros(1,FN);
            VLStep=(Max-Min)/(Te.VertLines);

            for f=1:FN
                FloorY(f)=Te.Surf(FloorX(f));
            end

            if Te.FloorLine==0 || ishandle(Te.FloorLine)==0
                % Draw horizontal line
                Te.FloorLine=line(FloorX,FloorY, 'LineWidth', 3*Te.LineWidth, 'Color', Te.FloorColor);

                % Draw vertical lines
                for v=1:Te.VertLines
                    Te.FloorVLx(v)=Min+v*VLStep;
                    Te.FloorVL(v)=line([Te.FloorVLx(v) Te.FloorVLx(v)-0.01],[Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-0.04],...
                                       'LineWidth', 2*Te.LineWidth, 'Color', Te.FloorColor);
                end
            else
                % Update horizontal line
                set(Te.FloorLine, 'XData', FloorX);
                set(Te.FloorLine, 'YData', FloorY);

                % Update vertical lines
                if Te.FloorVLx(1)<Min
                    Te.FloorVLx(1:end-1)=Te.FloorVLx(2:end);
                    Te.FloorVLx(end)=Te.FloorVLx(end-1)+VLStep;

                    for v=1:Te.VertLines
                        set(Te.FloorVL(v), 'XData', [Te.FloorVLx(v) Te.FloorVLx(v)-0.01]);
                        set(Te.FloorVL(v), 'YData', [Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-0.04]);
                    end
                end
                
                if Te.FloorVLx(end)>Max
                    Te.FloorVLx(2:end)=Te.FloorVLx(1:end-1);
                    Te.FloorVLx(1)=Te.FloorVLx(2)-VLStep;

                    for v=1:Te.VertLines
                        set(Te.FloorVL(v), 'XData', [Te.FloorVLx(v) Te.FloorVLx(v)-0.01]);
                        set(Te.FloorVL(v), 'YData', [Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-0.04]);
                    end
                end
            end        
        end
    end
end