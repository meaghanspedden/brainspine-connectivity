
ps=subject_m.pos;
figure
plot3(ps(:,1),ps(:,2),ps(:,3),'c.'); hold on
radius=bmodel.r;
origin=bmodel.o;

plot3(src_m.pos(:,1),src_m.pos(:,2),src_m.pos(:,3),'ro'); hold on
plot3(origin(1),origin(2),origin(3),'b*','markersize',12)

% Create a grid of points on the sphere's surface
theta = linspace(0, 2*pi, 100); % Azimuthal angle
phi = linspace(0, pi, 50);     % Polar angle
[theta, phi] = meshgrid(theta, phi);

% Calculate the Cartesian coordinates of the points on the sphere
x = origin(1) + radius * sin(phi) .* cos(theta);
y = origin(2) + radius * sin(phi) .* sin(theta);
z = origin(3) + radius * cos(phi);

% Create a 3D plot
surf(x, y, z,'facecolor','b','edgecolor','b','facealpha',0.1,'edgealpha',0.1);
axis equal; % Equal scaling for x, y, and z axes


