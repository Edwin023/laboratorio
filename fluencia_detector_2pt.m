function phi_2pt=fluencia_detector_2pt(sess_name,ndet,nx,ny,nz,det,dr,h)  


file_stats = dir([sess_name num2str(h) '_F.bin']);
    data_size = file_stats.bytes/(nx*ny*nz);
fid = fopen([sess_name num2str(h) '_F.bin'], 'rb');
    if(data_size==4)
        data_2pt = single(fread(fid, 'float32'));
    elseif(data_size==8)
        data_2pt = double(fread(fid, 'float64'));
    end
    fclose(fid);
    data_2pt = reshape(data_2pt(1:nx*ny), [nx,ny]);
    
    for i=1:ndet
       det(i,1)=det(i,1)/dr(1)+nx/2; 
       det(i,2)=det(i,2)/dr(2)+ny/2;
       
    end

    phi_2pt = zeros(ndet,1);
    for kk=1:ndet
    
        % Get the photon density from the detector voxel and neighboring voxels.
        % We assume detector radius of 1mm.
        i = round(det(kk,1));
        j = round(det(kk,2));
        

        phi_2pt(kk) = mean(abs([...
                                data_2pt(i-1, j-1) ...
                                data_2pt(i-1, j  ) ...
                                data_2pt(i-1, j+1) ...
                                data_2pt(i  , j-1) ...
                                data_2pt(i  , j  ) ...
                                data_2pt(i  , j+1) ...
                                data_2pt(i+1, j-1) ...
                                data_2pt(i+1, j  ) ...
                                data_2pt(i+1, j+1) ...
                               ]));
                         %  phi_2pt(kk) = abs(data_2pt(i,j ));
    end
    return