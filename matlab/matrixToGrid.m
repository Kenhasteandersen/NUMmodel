function [Cg,xg,yg,zg]=matrixToGrid(Cb,I,boxFile,gridFile,keepEmptyLayers,emptyLayerFillValue)

% function to rearrange boxes onto regular grid. No interpolation 
% is done. 
% [Cg,xg,yg,zg]=matrixToGrid(Cb,I,boxFile,gridFile,[keepEmptyLayers],[emptyLayerFillValue])
% rearranges vector of boxes Cb onto grid xg,yg,zg based on nominal box positions.
% I is index vector referenced to ALL boxes, so Cb=C_all_boxes(I); 
% e.g., if C is tracer at interior points (index vector Ii) use:
% [Cg,xg,yg,zg]=matrixToGrid(C,Ii,boxFile,gridFile,[keepEmptyLayers],[emptyLayerFillValue])
% The optional logical flag keepEmptyLayers (default is 0) determines 
% whether to eliminate empty layers before passing back the gridded field 
% (the default) or keep them. If keepEmptyLayers=1, returned empty layers 
% will simply be filled with the land (NaN)-sea (1) mask corresponding to the grid. 
% To use a value other than "1" for the mask pass the optional fill value emptyLayerFillValue. 
% To do 2d gridding: Suppose C is tracer at interior points but only 
% want to grid onto a single depth (zLayer):
% I=find(Zboxnom==zLayer);
% I1=find(Zboxnom(Ii)==zLayer); 
% [Cg,xg,yg,zg]=matrixToGrid(C(I1),I,boxFile,gridFile); 

if nargin<4
  error('ERROR: must pass 4 arguments');
end
if nargin<5
  keepEmptyLayers=0; % default is to DELETE empty layers
end  

load(gridFile,'x','y','z','nx','ny','nz')

if keepEmptyLayers
  load(gridFile,'bathy')
  if nargin<6 % specify default fill value for empty layer
    emptyLayerFillValue=1;
  end  
end

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','ixBox','iyBox','izBox','nb')  

nt=size(Cb,2);
xg=x;
yg=y;
if ~iscell(x) % regular GCM topology
  if isempty(I) % all boxes
	I=[1:nb]';
  end  
  Zb=Zboxnom(I);
  zg=unique(Zb);
  Cg=repmat(NaN,[nx ny nz nt]);
  tmp=repmat(NaN,[nx ny nz]);  
  idx=sub2ind([nx ny nz],ixBox(I),iyBox(I),izBox(I));
  for it=1:nt
	tmp(idx)=Cb(:,it);    
	Cg(:,:,:,it)=tmp;
  end  
% Only keep layers with any data
  if ~keepEmptyLayers
	kl=find(ismember(z,zg));
	Cg=Cg(:,:,kl,:);
  else
	kln=find(~ismember(z,zg));  
    bathy(bathy==0)=NaN;
    bathy(~isnan(bathy))=emptyLayerFillValue;
    for it=1:nt
	  Cg(:,:,kln,it)=bathy(:,:,kln);
    end
  end  
  Cg=squeeze(Cg);
else % assume LLC grid
  load(boxFile,'boxfacenum')  
  nFaces=length(x);
  if isempty(I) % all boxes
	I=[1:nb]';
  end
  Zb=Zboxnom(I);
  zg=unique(Zb);
  if ~keepEmptyLayers
	kl=find(ismember(z,zg));
  else
	kln=find(~ismember(z,zg));	
  end  
  Cg=cell(nFaces,1);
  for iF=1:nFaces	
	iif=find(boxfacenum(I)==iF); % indexed to I
	if ~isempty(iif)
	  Cg{iF}=repmat(NaN,[nx{iF} ny{iF} nz nt]);
	  tmp=repmat(NaN,[nx{iF} ny{iF} nz]);	  
	  Ibox=I(iif); % global box number
	  idx=sub2ind([nx{iF} ny{iF} nz],ixBox(Ibox),iyBox(Ibox),izBox(Ibox));
	  for it=1:nt
		tmp(idx)=Cb(iif,it);
		Cg{iF}(:,:,:,it)=tmp;
	  end
%     Only keep layers with any data
	  if ~keepEmptyLayers
		Cg{iF}=Cg{iF}(:,:,kl,:);
	  else
		bathy{iF}(bathy{iF}==0)=NaN;		
		bathy{iF}(~isnan(bathy{iF}))=emptyLayerFillValue;		
		for it=1:nt
		  Cg{iF}(:,:,kln,it)=bathy{iF}(:,:,kln);
		end
	  end
	  Cg{iF}=squeeze(Cg{iF});
	else
	  Cg{iF}=[];
	end  
  end  
end

