function varargout = ginput_ax(ha,n)
% function varargout = ginput_ax(ha,n)
%
% Alternative to ginput; restricts selection to specifed axis
%
% INPUTS    ha  axis handle
%           n   number of selection points, if empty defaults to 1
%
% Function copied from Matt Tearle, Mathworks, Aug 2021: 
% https://www.mathworks.com/matlabcentral/answers/7528-ginput-in-a-gui
%

if nargin<2
    n=1;
end
k = 0;
xy = zeros(n,2);
hf = get(ha,'parent');
figure(hf);
set(hf,'WindowButtonMotionFcn',@changepointer)
set(ha,'ButtonDownFcn',@getpoints)
hp = get(ha,'children');
ht = get(hp,'hittest');
set(hp,'hittest','off')
axlim = get(ha,'Position');
fglim = get(hf,'Position');
x1 = axlim(1)*fglim(3) + fglim(1);
x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
y1 = axlim(2)*fglim(4) + fglim(2);
y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
waitfor(hf,'WindowButtonMotionFcn',[])
if iscell(ht)
    for jj=1:length(ht)
        set(hp(jj),'hittest',ht{jj})
    end
else
    set(hp,'hittest',ht)
end
if nargout==2
    varargout{1} = xy(:,1);
    varargout{2} = xy(:,2);
else
    varargout{1} = xy;
end
      function changepointer(~,~)
          pntr = get(0,'PointerLocation');
          if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
              set(hf,'Pointer','crosshair')
          else
              set(hf,'Pointer','crosshair') %EK changed to crosshair
          end
      end
      function getpoints(hObj,~,~)
          cp = get(hObj,'CurrentPoint');
          k = k+1;
          xy(k,:) = cp(1,1:2);
          if k==n
              set(hf,'Pointer','arrow')
              set(hf,'WindowButtonMotionFcn',[])
              set(ha,'ButtonDownFcn',[])
          end
      end
end