function Cb=gridToMatrix(Cg,I,boxFile,gridFile,is2dArray)
% Transform grid variable Cg into vector of boxes represented by 
% (global) indices in I. Output is the same length as I. If I is empty, 
% the full 3-d grid variable is returned as a vector.

% Cg(x,y,z,[t]) or Cg(x,y,[t])
% If passing Cg(x,y,t) then is2dArray MUST be set to 1

if nargin<5
  is2dArray=0;
end

load(gridFile,'nx','ny','nz','gridType')
load(boxFile,'ixBox','iyBox','izBox','nb')  

if ~strcmp(gridType,'llc_v4') % regular GCM topology
% check
  if (size(Cg,1)~=nx) || (size(Cg,2)~=ny)
	error('ERROR: dimensions of Cg do not match expected size!')
  end

  if is2dArray % 2-d array with multiple time slices: interpret 3d dimension as time
    nt=size(Cg,3);
    nz1=1;
  else % interpret 3d dimension (if any) as depth and 4th dimension (if any) as time
	nz1=size(Cg,3);
	if (nz1~=nz) & (nz1~=1)
	  error('ERROR: can only pass a full 3-d grid or a single 2-d layer corresponding to boxes referenced by I')
	end
    nt=size(Cg,4);	
  end
  
  if isempty(I) % all boxes
	I=[1:nb]';
  end  

  if nz1==1 % a 2-d layer has been passed
    k=unique(izBox(I));
    if length(k)>1
      error('ERROR: Boxes referenced by I represent more than one layer when a single 2-d grid layer has been passed')
    end
  end
  
  if nz1==nz % a full 3-d grid has been passed
    idx=sub2ind([nx ny nz],ixBox(I),iyBox(I),izBox(I));
  else % a 2-d layer has been passed
    idx=sub2ind([nx ny],ixBox(I),iyBox(I));
  end
  
  Cb=zeros([length(idx) nt]);
  for it=1:nt
	if nz1==nz % a full 3-d grid has been passed
	  Cgtmp=Cg(:,:,:,it);
	else % a 2-d layer has been passed
	  Cgtmp=Cg(:,:,it);
	end  
	Cb(:,it)=Cgtmp(idx);
  end
  
%   if is2dArray
% 	for it=1:nt
% 	  Cgtmp=Cg(:,:,it);
% 	  Cb(:,it)=Cgtmp(idx);
% 	end    
%   else
% 	for it=1:nt
% 	  if nz1==nz % a full 3-d grid has been passed
% 		Cgtmp=Cg(:,:,:,it);
% 	  else % a 2-d layer has been passed
% 		Cgtmp=Cg(:,:,it);
% 	  end  
% 	  Cb(:,it)=Cgtmp(idx);
% 	end
%   end  
  
%   if nt==1 % only single time dimension passed
% 	Cb=Cg(idx);    
%   else % multiple time dimensions    
%     Cb=zeros([length(idx) nt]);
%     if is2dArray % 3d dimension interpreted as time
% 	  for it=1:nt
% 		Cgtmp=Cg(:,:,it);
% 		Cb(:,it)=Cgtmp(idx);
%       end    
%     else % 
% 	  for it=1:nt
% 		if nz1==nz % a full 3-d grid has been passed
% 		  Cgtmp=Cg(:,:,:,it);
% 		else % a 2-d layer has been passed
% 		  Cgtmp=Cg(:,:,it);
% 		end  
% 		Cb(:,it)=Cgtmp(idx);
% 	  end
% 	end  
%   end  
else % LLC grid
  nFaces=length(nx);
  
%   nFaces=length(Cg);

% check
  for iF=1:nFaces
	if (size(Cg{iF},1)~=nx{iF}) || (size(Cg{iF},2)~=ny{iF})
	  error('ERROR: dimensions of Cg do not match expected size!')
	end

	if is2dArray % 2-d array with multiple time slices: interpret 3d dimension as time
	  nt{iF}=size(Cg{iF},3);
	  nz1=1;
	else % interpret 3d dimension (if any) as depth and 4th dimension (if any) as time
	  nz1=size(Cg{iF},3);
	  if (nz1~=nz) & (nz1~=1)
		error('ERROR: can only pass a full 3-d grid or a single 2-d layer corresponding to boxes referenced by I')
	  end
	  nt{iF}=size(Cg{iF},4);
	end  
  end
  
  if any(cellfun(@(x)x==nt{1},nt)~=1)
    error('ERROR: All faces must have the same number of time slices')
  end
  nt=nt{1};
  
  load(boxFile,'boxfacenum')
  if isempty(I) % all boxes
	I=[1:nb]';
  end

  if nz1==1 % a 2-d layer has been passed
    k=unique(izBox(I));
    if length(k)>1
      error('ERROR: Boxes referenced by I represent more than one layer when a single 2-d grid layer has been passed')
    end
  end

  Cb=repmat(NaN,[nb nt]);
  for iF=1:nFaces	
	iif=find(boxfacenum(I)==iF); % indexed to I
	if ~isempty(iif)
	  Ibox=I(iif); % global box number
	  if nz1==nz % a full 3-d grid has been passed
		idx=sub2ind([nx{iF} ny{iF} nz],ixBox(Ibox),iyBox(Ibox),izBox(Ibox));
	  else % a 2-d layer has been passed
		idx=sub2ind([nx{iF} ny{iF}],ixBox(Ibox),iyBox(Ibox));
	  end
	  for it=1:nt
		if nz1==nz % a full 3-d grid has been passed
		  Cgtmp=Cg{iF}(:,:,:,it);
		else % a 2-d layer has been passed
		  Cgtmp=Cg{iF}(:,:,it);
		end  
		Cb(Ibox,it)=Cgtmp(idx);
	  end

% 	  if is2dArray
% 		for it=1:nt
% 		  Cgtmp=Cg{iF}(:,:,it);
% 		  Cb(Ibox,it)=Cgtmp(idx);
% 		end
% 	  else
% 		for it=1:nt
% 		  if nz1==nz % a full 3-d grid has been passed
% 			Cgtmp=Cg{iF}(:,:,:,it);
% 		  else % a 2-d layer has been passed
% 			Cgtmp=Cg{iF}(:,:,it);
% 		  end  
% 		  Cb(Ibox,it)=Cgtmp(idx);
% 		end
% 	  end  	  
    end
  end
  Cb=Cb(I,:);

%   Cb=repmat(NaN,[nb 1]);
%   for iF=1:nFaces	
% 	iif=find(boxfacenum(I)==iF); % indexed to I
% 	if ~isempty(iif)
% 	  Ibox=I(iif); % global box number
% 	  if nz1==nz % a full 3-d grid has been passed
% 		idx=sub2ind([nx{iF} ny{iF} nz],ixBoxGlob(Ibox),iyBoxGlob(Ibox),izBoxGlob(Ibox));
% 	  else % a 2-d layer has been passed
% 		idx=sub2ind([nx{iF} ny{iF}],ixBoxGlob(Ibox),iyBoxGlob(Ibox));
% 	  end
% 	  Cb(Ibox)=Cg{iF}(idx);
%     end
%   end
%   Cb=Cb(I);
end
