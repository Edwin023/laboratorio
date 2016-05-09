
close all

fid=fopen('sample1_F.bin','rb');
fluencia=fread(fid,'float32');
fluencia=reshape(fluencia,[240 240 153]);
%fluencia2D=reshape(log(fluencia(2:end,120,2:end)),[239 152])';
fluencia2D=reshape((fluencia(2:end,120,2:end)),[239 152])';
figure;
imagesc((fluencia2D).^(1/3))
colorbar;
xlabel('x(mm)','fontsize',16) 
ylabel('z(mm)','fontsize',16)


% %figure;
% subplot(1,2,2)
% imagesc(log(reshape(fluencia(2:end,50,2:end),[99 99])'))
% colorbar;
% xlabel('y(mm)','fontsize',16) 
% ylabel('z(mm)','fontsize',16)