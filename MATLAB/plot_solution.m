
% Generate Solution
hold on
[t, y] = ode45('lorenzSynthesized', [0 1], [0.55 .6 .6]);
cla

%[t, y] = ode45('lorenzHui', [0 500], [0.1 .5 1]);
%[t, y] = ode45('lorenzHui2', [0 500], [0.95 .5 .95]);
% Create 3d Plot
%[t, y] = ode45('lorentzPietro', [0 50], [0.5 .5 .8]);
%plot3(y(:,2), y(:,3), y(:,1))

zlabel('n1')
xlabel('n1c')
ylabel('n2c')
axis([0 1 0 1 0 1])
fixed=[]
%hold on
 opts = odeset('RelTol',1e-5,'AbsTol',1e-5)
for i=1:10
    for j=1:10  
        view(6*i, 40);   % Rotate the view
        
       % set(gca, 'XLim', [-20 20])
       % set(gca, 'YLim', [-25 25])
       % set(gca, 'ZLim', [0 50])
        
       %Remove the axes
        %axis off
        
        [t, y] = ode45('lorenzSynthesized', [0 10], [.5 (0+0.1*i) .1*j]);
        fixed=[fixed;y(end,:),y(end-20,:)]
        y(end,1)=nan;
        y(end,2)=nan;
        y(end,3)=nan;
        
        patch(y(:,2), y(:,3),y(:,1),y(:,1),'EdgeColor','interp','FaceColor','none')
       % Create 3d Plot
        %plot3(y(:,2), y(:,3), y(:,1),'marker','.');
        zlabel('n1')
        xlabel('n1c')
        ylabel('n2c')
        axis([0 1 0 1 0 1])
        hold on
%         holf 
%         Generate the file name for the new frame
%         name = sprintf('frames/frame-%02d.png', i);
%     
%     Save the current figure as a png file using the 
%     file name specified above. 
%     print('-dpng', name); 
    end
end

[t, y] = ode45('lorenzHui', [0 1], [(0+.2) (0+0.8) .8]);
y(end,1)=nan;
y(end,2)=nan;
y(end,3)=nan;

%patch(y(:,2), y(:,3),y(:,1),y(:,1),'EdgeColor','interp','FaceColor','none')

m=y;
[t, y] = ode45('lorenzHui', [0 1], [(0+.8) (0+0.2) .2]);
y(end,1)=nan;
y(end,2)=nan;
y(end,3)=nan;

%patch(y(:,2), y(:,3),y(:,1),y(:,1),'EdgeColor','interp','FaceColor','none')
ax = gca;
box on
ax.BoxStyle = 'full';

