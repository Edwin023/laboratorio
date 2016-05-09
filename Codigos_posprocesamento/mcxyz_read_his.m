function [phi,I]=mcxyz_read_his(sess_name, njobs)

% USAGE:
%
%    tMCimg_validate_output(sess_name, njobs)
%     
%    1. Change directory to where the tMCimg input and output files are stored.
%    2. Run tMCimg_validate_output giving the name of the tMCimg session and 
%       the number of jobs that the session was divided into. 
%
% DESCRIPTION:
%     
%    Validate the output of the tMCimg implementation of Monte Carlo. 
%
% EXAMPLE:
%
%    >> cd <tMCimg files root directory>
%    >> ls sample_*.2pt
%    sample_10.2pt  sample_2.2pt  sample_4.2pt  sample_6.2pt  sample_8.2pt
%    sample_1.2pt   sample_3.2pt  sample_5.2pt  sample_7.2pt  sample_9.2pt
%
%    >> tMCimg_validate_output('sample_', 10);
%
% AUTHORS:   Erin Buckley, (ebuckle2@gmail.com)
%            Jay Dubb,     (jdubb@nmr.mgh.harvard.edu)
%
   
%     global c; c=2.99792458e10;
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in the input data to tMCimg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    config = mcxyz_read_inp([sess_name num2str(1)]);
    %phot_num   = config.phot_num;
    %src        = config.srcpos;
   % time_gates = config.temporal_gates;
    %nx         = config.Dx; 
    %ny         = config.Dy; 
    %nz         = config.Dz;
    N_regions  = config.tiss_num;
    %mus        = config.tiss_prop(1,1); 
    %g          = config.tiss_prop(:,3);
    mua        = config.tiss_prop(:,1); 
    %det        = config.detpos;
    ndet       = config.num_det;

%     %musp = (1-g)*mus;
%     t0 = time_gates(1);
%     tt = time_gates(2);
%     delta_t = time_gates(3);

    % Calculate nt (number of time gates) taking into account floating point 
    % division errors
%     nt_float = (tt-t0)/delta_t;
%     nt_int   = double(int8((tt-t0)/delta_t));
%     delta_t_r = abs(nt_float - nt_int) * delta_t;
%     delta_t_too_small = 1e-6 * delta_t;
%     if(delta_t_r < delta_t_too_small)
%         nt = nt_int;
%     else
%         nt = ceil(nt_float);
%     end 

    %r = sqrt((src(1,1)-det(1,1))^2 + (src(1,2)-det(1,2))^2 + (src(1,3)-det(1,3))^2);
%     v = c/n;

    % Get the surface plane of the homogeneous medium; we assume the detectors are 
    % all on the homogeneous plane.
    %z = det(1,3);

    % Get tissue absorption and replace tissue voxels in seg file 
    % with the tissue's absortion value
%     slab_seg = load_bin(config.seg_file);
%     mua_vol = single(zeros([size(slab_seg) nt_int]));
%     i_slab_seg = find(slab_seg == 1);
%     mua_vol(i_slab_seg) = mua;
%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get fluence from history data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi = zeros(njobs,ndet);
    for h=1:njobs

        fid = fopen([sess_name num2str(h) '.his'], 'rb');
        data_his = fread(fid, inf, 'float32');
        nhits = size(data_his, 1) / (2*N_regions+1);
        data_his = reshape(data_his, [2*N_regions+1 nhits])';
        fclose(fid);

        for k=1:ndet
            % Get path lengths of all photons that hit the kth detector
            %L_k = data_his(find(data_his(:,1) == k-1), 2:end);
            L_k = data_his(find(data_his(:,1) == k-1), 2:1+N_regions);
            
            % Get the total time in seconds each photon spent traveling before 
            % reaching the kth detector
%             T_k = sum(L_k, 2)/v;


            % Calculate fluence over time t, for detector k; implementation of
            % formula for obtaining fluence over time from history data:
            % 
            %           Sigma{ i=0:N_photons(t) } ( Pi{ j=1:N_regions } (exp(-mua(j)*L(i,j))) )
            % phi(t) =  -------------------------------------------------------------------
            %                               N_photons(t)
            %
            % Reference: "Three dimensional Monte Carlo code for photon migration 
            % through complex heterogeneous media including the adult human head"
            % D. A. Boas, J. P. Culver, J. J. Scott, A. K. Dunn 
            %
           % for t=1:nt
                % Get all photons hitting det k at time gate t_i
            %    t_i = (t-1)*delta_t;
             %   ii = find(T_k > t_i & T_k <= (t_i+delta_t));

                % L: path lengths of all photons that hit the kth detector at time t
               % L = L_k(ii,:);
                L = L_k;

                % Calculate fluence at detector k, that is the sum of 
                % the fluences of all the photons that hit detector k
                N_photons = size(L, 1);
                for i=1:N_photons

                    % Calculate fluence for 1 photon 
                    w = 1;
                    for j=1:N_regions
                        w = w * exp(-mua(j)*L(i,j));
                    end
                    phi(h, k) = phi(h, k) + w;
                end
                phi(h, k) = phi(h, k) ;
                I(h, k) = N_photons;
            %end
        end
    end
    clear data_his

%     phimean_his  = squeeze(mean(phi,1));
%     phistd_his   = squeeze(std(phi,0,1));
%     Imeans       = squeeze(mean(I,1));
%     Istd         = squeeze(std(I,0,1));