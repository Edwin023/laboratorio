function mcxyz_ejecutor(mc_sess, time_min,time_prop, dims, optpos, tiss_prop, num_jobs)

%
% USAGE:
%
%    tMCimg_validate(mc_sess, num_phot, time_gates, dims, optpos, tiss_prop, num_jobs, platform)
%
% DESCRIPTION:
%
%    This function reads the input files for a monte carlo simulation using tMCimg, runs tMCimg and
%    then validates the output using the tMCimg_validate_output.m function. 
%
%
% INPUTS:
%
%    mc_sess    - Name of the tMCimg session; this determines the names of all the 
%                 .inp, .2pt and .his files 
%    
%    num_phot   - The number of photons to launch. 
%    time_gates - Vector containing the start and end times during which to record the 
%                 photon fluence and the time step size into which to break up the time gate.  
%    dims       - Dimensions of the volume containging the homogeneous medium
%    optpos     - Positions of one light source and detectors. The source is always the first entry. 
%    tiss_prop  - A struct array containing the following tissue properties fields
%        scattering
%        absorption
%        refraction
%        anisotropy
%    num_jobs    - The number of concurrent jobs to divide the number of photons into for a cpu cluster. 
%    platform    - The operating system you're running this script on. The choices are: 
%                  'windows', 'linux', or 'mac'. If this argument is not provided, then 
%                  the function will pause and ask the user to run tMCimg manually, then to press
%                  'Enter' when the output has been generated.  
%
% EXAMPLE:
%
%    The following example generates a 100x100x100 homogeneous, semi-infinite volume, places 
%    5 optodes on the surface - one source and 5 detectors all spaced 10mm away from the source. 
%    Then in the last step, the output of the tMCimg monte carlo simulation is validated.
%
%
%    cd /<root path>/montecarlo/example/validation
%    tiss_prop(1).scattering = 50;
%    tiss_prop(1).absorption = 0.0100;
%    tiss_prop(1).refraction = 1;
%    tiss_prop(1).anisotropy = 0.9800;
%    optpos = [50 50 10;  50 40 10;   50 60 10;  40 50 10;  60 50 10];
%    tMCimg_validate('sample', 10000000, [0 1.0e-9 0.1e-9], [100 100 100], optpos, tiss_prop, 10);
%
%
% AUTHORS:  Jay Dubb,     (jdubb@nmr.mgh.harvard.edu)
%           Erin Buckley, (ebuckle2@gmail.com)
%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check arguments if they exists - if not supply default values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(~exist('mc_sess') | isempty(mc_sess))
        mc_sess='sample';
    end

    if(~exist('time_min') | isempty(time_min))
        num_phot=10000000;
    end

  

    if(~exist('dims') | isempty(dims))
        dims = [100 100 100];
    end

    

    % Dimensions
    Dx = dims(1);
    Dy = dims(2);
    Dz = dims(3);

    if(~exist('optpos') | isempty(optpos))

        % Default radius is 10mm
        r=10;

        % Default spacing of detector around source at 
        % 45 degree intervals - i.e., 8 detectors
        rho=45;

        % Source positions
        src(1,:) = [0   0  0.1];    
        dets=find_detectors(src(1,:), r, rho);
        optpos=[src; dets];

    end
    a = optpos(1,3)/0.1;

    if(~exist('tiss_prop'))
        tiss_prop = set_tiss([], 'other', 50, .01, 1, .98, 1);
    end

    if(~exist('num_jobs') | isempty(num_jobs))
        num_jobs = 10; 
    end
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate semi-infinite homogeneous volume 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slab_seg = zeros(Dx, Dy, Dz);
    slab_seg(:,:,(a+1):end) = 1;
    %slab_seg(:,:,(a+9):end) = 2;
    slab_seg=reshape(slab_seg,[Dx*Dy*Dz 1]);

    % Create .bin file from segemented head
    seg_fn = strcat(mc_sess, [num2str(num_jobs) '_T.bin']);
    fid = fopen(seg_fn, 'wb');
    fwrite(fid, slab_seg, 'uint8');
    fclose(fid);

    % % Save text file with the dimensions
    %i = strfind(seg_fn, '.');
    %seg_fn_dims = [seg_fn(1:i(end)-1) '_dims.txt'];
   % dims = size(slab_seg); 
    %save(seg_fn_dims, 'dims', '-ascii');
    
    tam_bin=[0.1 0.1 0.1];
    mcflag=2;
    boundary=1;
    rfocus=[0 0 inf];
    par_source=[0.2 0.01];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate input files for tMCimg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     gen_mcxyz_valid_input(mc_sess, ...
                               dims, ...
                                tam_bin,...
                                mcflag,...
                                boundary,...
                                optpos, ...
                                rfocus,...
                                par_source,...
                               time_min, ...
                               time_prop,...
                               tiss_prop, ...
                                num_jobs)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run tMCimg from matlab if the platform argument is supplied;
    % otherwise just pause and let user run it manually, then 
    % have user tell us that tMCimg is done by pressing 'Enter'.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    display(['Generated input files in ' pwd]);
%     if(isempty(platform))
%         display('Now run the Monte Carlo forward simulation using tMCimg. (The csh script tMCimg_all.csh  '); 
%         display('(The csh scripts tMCimg_all.csh and tMCimg_all_pbsubmit.csh run the tMCimg jobs for all optodes.)');
%         c = input('When it''s done, press ENTER to validate tMCimg output or ''q'' to quit...: ', 's');
%         if(c == 'q')
%             return;
%         end
%     else
%         if(strcmp(platform, 'linux'))
%             system('chmod 755 ./tMCimg_all.sh');
%             system('./tMCimg_all.sh');
%         elseif(strcmp(platform, 'windows'))
%             system('tMCimg_all.bat');
%         end
%     end   
        system('mcxyz_all.bat');