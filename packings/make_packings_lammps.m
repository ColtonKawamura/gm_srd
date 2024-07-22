% The purpose of this script is the create data files
% for LAMMPS granular packing simulations. Using MATLAB
% or OCTAVE, you can call this function using 
% MakeLammpsICs(N,W) in the command line where you
% determine N and W as numbers. 
%
% W = width of square particle top and bottom boundaries
% defined by number of particles EXCLUDING the particles
% at boundary. For example, a W = 5 will create a 4x4 wall.
%
% N = number of particles that are between the walls.
%
% Example, a "MakeLammpsICs(1000,5)" call will create a 
% packing with 4x4x2=32 particles on the top and bottom walls
% with 1000 particles between them for a total of 1032 
% particles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function make_packings_lammps(N,W)
%% Set Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%  N=50; W=5;

varsname=['N' num2str(N) 'W' num2str(W)];

% Names of files
filename=['./data_IC' varsname '.granular']; 

% Average Particle size in Simulation Units
d_ave=1;

% Density of particles
rho=1;



%% Calculate Values Needed for Text File

% Size and density of wall particles
    rho_wall=0.05*rho;
    d_max=1.3*d_ave;
    d_wall=2.6*d_ave; %CAK Changed from 1.6
    
% Calculate Width of Box in Simulation Units
    Wbox=ceil(W.*d_ave);
    

% Wall particle position vectors
% All wall particles are d_max
    sq_wall_pos=[d_max/2:d_max:ceil(W*d_ave)];

    all_x_wall=ones(length(sq_wall_pos),1)*sq_wall_pos;
    all_x_wall=all_x_wall(:)'-(W*d_ave)/2;
%     all_x_wall=[all_x_wall,all_x_wall+0.3*d_max,all_x_wall-0.3*d_max];

    all_y_wall=(ones(length(sq_wall_pos),1)*sq_wall_pos)';
    all_y_wall=all_y_wall(:)'-(W*d_ave)/2;
%     all_y_wall=[all_y_wall,all_y_wall+0.3*d_max,all_y_wall-0.3*d_max];

    
% Calculate number of particles in tightly packed walls
% Nlower=ceil((Wbox/d_max).^2);
% Nupper=ceil((Wbox/d_max).^2);
    Nlower=length(all_y_wall);
    Nupper=Nlower;

    Ntot=N+Nlower+Nupper;

    
% Flow particles per layer
    sq_flow=[(1.5*d_ave):(1.5*d_ave):Wbox];
    
    
%Placement of flow particles
    layer_x_flow=ones(length(sq_flow),1)*sq_flow;
    layer_x_flow=layer_x_flow(:)'-(W*d_ave)/2;

    layer_y_flow=(ones(length(sq_flow),1)*sq_flow)';
    layer_y_flow=layer_y_flow(:)'-(W*d_ave)/2;
    
    
% Approxmate number of flow layers needed
    FlowLayers=ceil(N/length(layer_y_flow));

    
% Set Height of Box
    %Hbox=(10*d_max)+ceil(N./((W/3.2).^2))*3;
    Hbox=(7*d_max)+(3*(FlowLayers+2)*d_ave);

    
% Raise and lower all particles within 0.5 grain diameters to create
% roughness on walls
    add_lower=0.2*randn(1,Nlower);
    add_upper=0.2*randn(1,Nupper);

    z_wall_upper_pos=Hbox-2.5*d_max+add_upper;%[Hbox-2.5*d_max+add_upper,Hbox-d_max+zeros(1,(Nupper/3)),Hbox-1.5*d_max+zeros(1,(Nupper/3))];
    z_wall_lower_pos=add_lower+1*d_max;%[add_lower+1*d_max,zeros(1,Nlower/3),zeros(1,Nlower/3)-0.5*d_max];
    %z_wall_upper_pos=[Hbox-2.5*d_max+add_upper,Hbox-1.5*d_max+add_upper,Hbox-0.5*d_max+add_upper];


% Set Positions for all flow particles
    sq_flow_h=[(2*d_max):(1*d_ave):(Hbox-2*d_max)];
    
    all_z_flow=ones(length(layer_x_flow),1)*sq_flow_h;
    all_z_flow=all_z_flow(:)';
    all_z_flow = all_z_flow + 0*(rand(size(all_z_flow))-0.5);
    
    all_x_flow=(ones(length(sq_flow_h),1)*layer_x_flow)';
    all_x_flow=all_x_flow(:)';
    all_x_flow = all_x_flow + 0.5*(rand(size(all_z_flow))-0.5);

    all_y_flow=(ones(length(sq_flow_h),1)*layer_y_flow)';
    all_y_flow=all_y_flow(:)';
    all_y_flow = all_y_flow + 0.5*(rand(size(all_z_flow))-0.5);
    
    % Move wall down
    
    z_wall_upper_pos = z_wall_upper_pos + 2 - (min(z_wall_upper_pos)-max(all_z_flow(1:N)));

%% Make Text File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID=fopen(filename,'w');

fprintf(fileID,['LAMMPS 3d granular data file\n'...
    '\n'...
    num2str(Ntot) ' atoms\n'...
    '0 bonds\n'...
    '0 angles\n'...
    '0 dihedrals\n'...
    '0 impropers\n'...
    '\n'...
    '3 atom types\n'...
    '0 bond types\n'...
    '0 angle types\n'...
    '0 dihedral types\n'...
    '0 improper types\n'...
    '\n'...
    num2str(-Wbox/2) ' ' num2str(Wbox/2) ' xlo xhi\n'...
    num2str(-Wbox/2) ' ' num2str(Wbox/2) ' ylo yhi\n'...
    num2str(-2*d_ave) ' ' num2str(Hbox) ' zlo zhi\n'...
    '\n'...
    '\n'...
    'Atoms\n'...
    '\n']);


atomID=1:Ntot;
atomtype=[ones(1,Nlower)*2,ones(1,Nupper)*3,ones(1,N)];
diameter=[ones(1,Nlower)*d_wall,ones(1,Nupper)*d_wall,0.2*randn(1,N)+d_ave];
rho_all=[ones(1,Nlower)*rho_wall,ones(1,Nupper)*rho_wall,ones(1,N)*rho];
x_all=[all_x_wall,all_x_wall,all_x_flow(1:N)];
y_all=[all_y_wall,all_y_wall,all_y_flow(1:N)];
z_all=[z_wall_lower_pos,z_wall_upper_pos,all_z_flow(1:N)];

Amat = [atomID;atomtype;diameter;rho_all;x_all;y_all;z_all];

%Format is atom_ID, atom_type, diameter, rho, x, y, z
formatSpec = '%.0f %.0f %.3f  %.3f %.4f %.4f %.4f\n';
fprintf(fileID,formatSpec,Amat);

fclose(fileID);

figure,
plot3(x_all,y_all,z_all,'.')
computational_time = toc;
fprintf('Elapsed time: %.2f seconds\n', computational_time)
end