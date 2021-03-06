% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate figure(1) illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
% Steven L. Jacques, March 3, 2013.

clear
format compact
clc
home
%global tissue

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci 
                        % 0 = don't save. Just check the program.

myname      = 'example1';% name for files: myname_T.bin, myname_H.mci  
nm          = 532;   	% desired wavelength of simulation
time_min    = 5;      	% time duration of the simulation
Nbins       = 100;    	% # of bins in each dimension of cube 
binsize     = 0.1; 	% size of each bin, eg. [cm] or [mm]

% Set Monte Carlo launch flags
mcflag      = 2;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 

% boundary condition
boundary    = 0;        % 0 = finite boundaries. Photons escape if outside cube.
                        % 1 = infinite medium. Let photons wander outside cube.
                        
% Sets position of launch
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.1;     % z of source

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 (uniform or Gaussian beam, not isotropic pt.)
radius      = 0.2;      % 1/e radius of beam at tissue surface
waist       = 0.010;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.7;      % trajectory projected onto x axis
uy0         = 0.4;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n\n-------- myname = %s ----------\n',myname))

if boundary==0
    disp('boundary is finite. Photon escapes outside tissue structure.')
else
    disp('boundary is infinite. Photon allowed to wander outside tissue structure.')
end
disp(' ')

%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

Nt=1;
for i=1:Nt
muav(i)=0.05;
musv(i)=1;
gv(i)=0;
end

Nd=2; % # of detectors
radDet=0.1; %radius od detectors

posDet=[0.6 0.5 0; 0.7 0.5 0  ];


% % Create tissue properties
% nm = 532; % wavelength
% tissue = makeTissueList(nm); % also --> global tissue(1:Nt).s
% Nt = length(tissue);
% for i=1:Nt
%     muav(i)  = tissue(i).mua;
%     musv(i)  = tissue(i).mus;
%     gv(i)    = tissue(i).g;
% end
% 
% Specify Monte Carlo parameters    
Nx = Nbins;
Ny = Nbins;
Nz = Nbins;
dx = binsize;
dy = binsize;
dz = binsize;
%% scale axes
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

% if isinf(zfocus), zfocus = 1e12; end
% 
% % %%%%%%
% % % CREATE TISSUE STRUCTURE T(y,x,z)
% % %   Create T(y,x,z) by specifying a tissue type (an integer)
% % %   for each voxel in T.
% % %
% % %   Note: one need not use every tissue type in the tissue list.
% % %   The tissue list is a library of possible tissue types.
% % 
T = double(zeros(Ny,Nx,Nz)); 
% % 
% % T = T + 4;      % fill background with skin (dermis)
% % 
 zsurf = 0.100;  % position of air/skin surface
% % 
for iz=1:Nz % for every depth z(iz)
% % 
% %     % air
     if iz>round(zsurf/dz)
         T(:,:,iz) =1; 
     end
% % 
% %     % epidermis (60 um thick)
% %     if iz>round(zsurf/dz) & iz<=round((zsurf+0.0060)/dz)
% %         T(:,:,iz) = 5; 
% %     end
% % 
% %     % blood vessel @ xc, zc, radius, oriented along y axis
% %     xc      = 0;            % [cm], center of blood vessel
% %     zc      = Nz/2*dz;     	% [cm], center of blood vessel
% %     vesselradius  = 0.0500;      	% blood vessel radius [cm]
% %     for ix=1:Nx
% %             xd = x(ix) - xc;	% vessel, x distance from vessel center
% %             zd = z(iz) - zc;   	% vessel, z distance from vessel center                
% %             r  = sqrt(xd^2 + zd^2);	% r from vessel center
% %             if (r<=vesselradius)     	% if r is within vessel
% %                 T(:,ix,iz) = 3; % blood
% %             end
% % 
% %     end %ix
% %     
end % iz
% % 
% % 
% % %%
%  if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,boundary);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
        end
        fprintf(fid,'%d\n',Nd);
        fprintf(fid,'%0.4f\n',radDet);
        for i=1:Nd
           for j=1:3
               fprintf(fid,'%0.4f\n',posDet(i,j));
           end
        end
    fclose(fid);

    %% write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

     toc
 %end % SAVEON
% 
% %% Look at structure of Tzx at iy=Ny/2
% Txzy = shiftdim(T,1);   % Tyxz --> Txzy
% Tzx  = Txzy(:,:,Ny/2)'; % Tzx
% 
% %% Tissue Structure
% figure(1); clf
% fsz = 18;  % font size 
% imagesc(x,z,Tzx,[1 Nt])
% hold on
% set(gca,'fontsize',fsz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% colorbar
% cmap = mycolormap2(Nt);
% colormap(cmap)
% set(colorbar,'fontsize',1)
% % label colorbar
% for i=1:Nt
%     yy = zmin + (Nt-i)/(Nt-1)*zdiff;
%     text(xmin + xdiff*1.13,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',12)
% end
% text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)
% axis equal image
% axis([xmin xmax zmin zmax])
% 
% 
% % draw launch
% N = 10; % # of beam rays drawn
% switch mcflag
%     case 0 % uniform
%         for i=0:5
%             for j=-2:2
%             plot( [xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
%             plot(-[xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
%             end
%         end
% 
%     case 1 % Gaussian
%         for i=0:5
%             for j=-2:2
%             plot( [xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
%             plot(-[xs+radius*i/5 xfocus + waist*j/2],[zs zfocus],'r-')
%             end
%         end
% 
%     case 2 % iso-point
%         for i=1:20
%             th = (i-1)/19*2*pi;
%             xx = Nx/2*cos(th) + xs;
%             zz = Nx/2*sin(th) + zs;
%             plot([xs xx],[zs zz],'r-')
%         end
% end

disp('done creating tissue.')

