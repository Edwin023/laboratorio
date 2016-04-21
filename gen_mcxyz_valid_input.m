function gen_mcxyz_valid_input(mc_sess, ...
                               seg_dims, ...
                                tam_bin,...
                                mcflag,...
                                boundary,...
                                optpos, ...
                                rfocus,...
                                par_source,...
                               Nphotons, ...
                               tiss_prop)

% USAGE: 
%
% gen_tMCimg_valid_input(mc_sess, ...
%                        seg_file, ...
%                        seg_dims, ...
%                        optpos, ...
%                        num_phot, ...
%                        time_gates, ...
%                        tiss_prop, ...
%                        num_jobs)
%
% DESCRIPTION:
%
% Function creates the .inp files for validating tMCimg simulation of light diffusion 
% in semi-infinite, homogeneous media. 
%
% EXAMPLE:
% seg = zeros(100,100,100);
% seg(:,:,10:end) = 1;
% save_bin('seg.bin', seg)
% tiss_prop(1).scattering = 50;
% tiss_prop(1).absorption = 0.0100;
% tiss_prop(1).refraction = 1;
% tiss_prop(1).anisotropy = 0.9800;
% optpos = [50 50 10;  50 40 10;   50 60 10;  40 50 10;  60 50 10];
% gen_tMCimg_valid_input('sample', ...
%                        './seg.bin', ...
%                        [100 100 100], ...
%                        optpos, ...
%                        1000000, ...
%                        [0 1.0e-9 0.1e-9], ...
%                        tiss_prop, ...
%                        10, ...); 
%
% Created Feb 20, 2008 by Jay Dubb, (jdubb@nmr.mgh.harvard.edu)
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Retrieve argument values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Nx = seg_dims(1);
    Ny = seg_dims(2);
    Nz = seg_dims(3);

    binsizex     = tam_bin(1); 	% size of each bin, eg. [cm] or [mm]
    binsizey     = tam_bin(2); 	% size of each bin, eg. [cm] or [mm]
    binsizez     = tam_bin(3); 	% size of each bin, eg. [cm] or [mm]binsizex     = tam_bin1(); 	% size of each bin, eg. [cm] or [mm]

% % Set Monte Carlo launch flags
% mcflag      = 2;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 

% boundary condition
%boundary    = 0;        % 0 = finite boundaries. Photons escape if outside cube.
                        % 1 = infinite medium. Let photons wander outside cube.
                        
% Sets position of launch
xs          = optpos(1,1);      	% x of source
ys          =optpos(1,2);       % y of source
zs          = optpos(1,3);     % z of source

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = rfocus(1);        % set x,position of focus
yfocus      = rfocus(2);        % set y,position of focus
zfocus      = rfocus(3);     	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 (uniform or Gaussian beam, not isotropic pt.)
radius      = par_source(1);      % 1/e radius of beam at tissue surface
waist       = par_source(2);  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0         = 0.0;      % trajectory projected onto x axis
uy0         = 0.0;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1
%%%%%%%%%%%%%%%%%%%%%%%%%

% % disp(sprintf('\n\n-------- myname = %s ----------\n',myname))
% % 
% % if boundary==0
% %     disp('boundary is finite. Photon escapes outside tissue structure.')
% % else
% %     disp('boundary is infinite. Photon allowed to wander outside tissue structure.')
% % end
% % disp(' ')

%%%%%%%%%% 
% Prepare Monte Carlo 
%%%



radDet=0.1; %radius od detectors

posDet=optpos(2:end,:);
posDet(1:end,1)=posDet(1:end,1);
posDet(1:end,2)=posDet(1:end,2);
Nd = size(posDet, 1);
 Nt = length(tiss_prop);
  t = tiss_prop;
  
    fid_mcxyz_master_win = fopen('mcxyz_all.bat', 'wt');
   % mci_files = cell(1, 1);
    %for j=1:num_jobs
        %mci_file = [mc_sess num2str(j)];
        mci_file = mc_sess;

        %fprintf(fid_tMCimg_master_lin, 'time tMCimg %s\n', inp_file); 
        %fprintf(fid_tMCimg_master_win, 'tMCimg.exe %s\n', inp_file); 
        fprintf(fid_mcxyz_master_win, 'mcxyz.exe %s\n', mci_file);

        %%%% Append extension to file name to create file. 
       % mci_files = {[mci_file '_H.mci']};
    %end
   % for j=1:num_jobs
        fid = fopen([mci_file '_H.mci'], 'wt');
   
     % run parameters
        fprintf(fid,'%d\n',Nphotons);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n', binsizex);
        fprintf(fid,'%0.4f\n', binsizey);
        fprintf(fid,'%0.4f\n', binsizez);
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
            fprintf(fid,'%0.9f\n',t(i).muav);
            fprintf(fid,'%0.9f\n',t(i).musv);
            fprintf(fid,'%0.9f\n',t(i).gv);
        end
        fprintf(fid,'%d\n',Nd);
        fprintf(fid,'%0.4f\n',radDet);
        for i=1:Nd
           for j=1:3
               fprintf(fid,'%0.4f\n',posDet(i,j));
           end
        end
    
   % end
   fclose(fid);
    
    %%Tmcimg
 
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Generate file names by appending the optode name and 
    %%%% .inp extension and then open the file for IO.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Create and write the data to each .inp file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    
    fclose(fid_mcxyz_master_win);
