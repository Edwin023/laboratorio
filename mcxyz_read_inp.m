function config =mcxyz_read_inp(mc_sess)

% USAGE: 
%
%     config = tMCimg_read_inp(mc_sess)
%
% DESCRIPTION:
%
%     Reads in tMCimg configuration data from the first .inp input file
%     into the structure config. The structure fields are
%
%         phot_num   
%         seed
%         srcpos
%         srcdirection
%         temporal_gates
%         seg_file
%         Dx
%         Dy
%         Dz
%         Sx_min
%         Sx_max
%         Sy_min
%         Sy_max
%         Sz_min
%         Sz_max
%         tiss_num
%         tiss_prop
%         num_det
%         det_rad
%         detpos
%
%     The input file name that it looks for is <mc_sess>.s1.inp, where mc_sess 
%     is the input argument. 
%
% EXAMPLE 1: 
%
%     >> dir slab_phi_0.*.inp
% 
%     slab_phi_0.d1.inp  slab_phi_0.d3.inp  slab_phi_0.d5.inp  slab_phi_0.d7.inp  slab_phi_0.d9.inp  slab_phi_0.s2.inp  slab_phi_0.s4.inp  
%     slab_phi_0.d2.inp  slab_phi_0.d4.inp  slab_phi_0.d6.inp  slab_phi_0.d8.inp  slab_phi_0.s1.inp  slab_phi_0.s3.inp  
%
%     >> config = tMCimg_read_inp('slab_phi_0');
%     >> config
% 
%     config = 
% 
%               phot_num: 10000000
%                   seed: -39111473
%                 srcpos: [40 60 10]
%           srcdirection: [0.2300 -0.1900 0.9500]
%         temporal_gates: [0 5.0000e-09 5.0000e-09]
%               seg_file: 'slab_seg.bin'
%                     Dx: 100
%                     Dy: 100
%                     Dz: 100
%                 Sx_min: 1
%                 Sx_max: 100
%                 Sy_min: 1
%                 Sy_max: 100
%                 Sz_min: 1
%                 Sz_max: 100
%               tiss_num: 1
%              tiss_prop: [0.6600 1.0000e-03 0.0190 1]
%                num_det: 12
%                det_rad: 1
%                 detpos: [12x3 double]
%        
%
% AUTHOR: Jay Dubb (jdubb@nmr.mgh.harvard.edu)
% DATE:   06/24/2009


    % Define and allocate structure for the output
    config = struct('time_min',0,'Dims',[0 0 0],'bins', [0 0 0] , 'mcflag', 0, 'boundary', 0,...,
    'srcpos',[0 0 0], 'focus',[0 0 0] , 'srcdirection',[0 0 0], 'par_sourse',[0 0], ...,
                      'tiss_num',0, 'tiss_prop',[],'num_det',0,'det_rad',0, 'detpos', []);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in parameters from the inp file for laser 1 
    % (that's always garanteed to be there). 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(mc_sess, 'rb');
    if(fid == -1)
        mc_sess_inp_file = strcat(mc_sess, '_H.mci');
        fid = fopen(mc_sess_inp_file, 'rb');
        if(fid == -1)
            mc_sess_inp_file = strcat(mc_sess, '.s1_H.mci');
            fid = fopen(mc_sess_inp_file, 'rb');
            if(fid == -1)
                config=[];
            end
        end
    end


    config.time_min = str2num(fgetl(fid));
    for i=1:3
    config.Dims(1,i) = str2num(fgetl(fid));
    end
    for i=1:3
    config.bins(1,i) = str2num(fgetl(fid));
    end
    
    
    config.mcflag=str2num(fgetl(fid));
    config.boundary=str2num(fgetl(fid));
    for i=1:3
    config.srcpos(1, i) = str2num(fgetl(fid));
    end
%     config.time_pop(1) = str2num(fgetl(fid));
%     config.time_pop(2) = str2num(fgetl(fid));
%     config.time_pop(3) = str2num(fgetl(fid));
%     
    for i=1:3
    config.focus(1, i) = str2num(fgetl(fid));
    end
    for i=1:3
    config.srcdirection(1,i) = str2num(fgetl(fid));
    end
    for i=1:2
    config.par_sourse(1,i) = str2num(fgetl(fid));
    end
    % Get tissue properties
    config.tiss_num = str2num(fgetl(fid));
    for i=1:config.tiss_num
        for j=1:3
        config.tiss_prop(i,j) = str2num(fgetl(fid));
        end
    end
    
    % Get rest of the optode positions
   config.num_det= str2num(fgetl(fid));
   
    
    if(config.num_det ~= 0)
       config.det_rad= str2num(fgetl(fid));
        for i=1:config.num_det
            for j=1:3
            config.detpos(i, j) = str2num(fgetl(fid));
            end
        end
    end
    fclose(fid);


