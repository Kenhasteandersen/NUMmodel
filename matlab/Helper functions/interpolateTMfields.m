%
% Interpolate the initial nutrient fields (N0, Si0) and the light field (parday)
% from one set of transport matrices onto the grid of another one.
%
% This is for transport matrices that come without those fields, such as
% UVicOSUpicdefault. The proper source for N and Si is the world ocean atlas
% through importN, but that needs the WOA netcdf files; this is what to do when
% they are not at hand.
%
% Three things need care, and are the reason this is not a plain interp3:
%
%  - The source fields have NaNs over land (the 2.8 degree parday is 46 % NaN).
%    Interpolating them directly spreads the NaNs into the ocean, so the gaps are
%    first filled from the nearest cell that has data. This is what the mask in
%    create_parday does. Note that zeros are kept: a zero in parday is the polar
%    night, not missing data, and masking zeros out fills the dark part of the
%    year with light from further south.
%
%  - The target grid reaches further towards the poles and deeper than the source
%    does, so the interpolation is linear inside the source grid and nearest
%    neighbour outside it.
%
%  - parday is a time series with one column per transport step, so it also has
%    to be interpolated in time, from 365/dtTransport columns of the source to
%    however many the target needs. The interpolation wraps around new year.
%
% In:
%  nTMtarget - which transport matrices to write the fields for, see
%              parametersGlobal
%  nTMsource - which to interpolate from (default 2, MITgcm_ECCO)
%  options.bOverwrite - whether to overwrite files that are already there
%              (default false)
%  options.fields - which fields to do (default all three)
%
% Example, to give the UVic matrices their fields:
%  interpolateTMfields(3)
%
function interpolateTMfields(nTMtarget, nTMsource, options)

arguments
    nTMtarget {mustBeInteger}
    nTMsource {mustBeInteger} = 2; % MITgcm_ECCO
    options.bOverwrite logical = false;
    options.fields cell = {'N0','Si0','parday'};
end

if nTMtarget == nTMsource
    error('The source and the target are the same transport matrices.')
end
%
% A setup is only needed to be able to call parametersGlobal, which is where the
% paths to the fields live:
%
pDummy = setupGeneralistsSimpleOnly;
pS = parametersGlobal(pDummy, nTMsource);
pT = parametersGlobal(pDummy, nTMtarget);
fprintf('Interpolating from %s onto %s\n', pS.TMname, pT.TMname);

gS = load(pS.pathGrid, 'x','y','z');
gT = load(pT.pathGrid, 'x','y','z');
fprintf('  source grid %dx%dx%d, target grid %dx%dx%d\n', ...
    numel(gS.x), numel(gS.y), numel(gS.z), numel(gT.x), numel(gT.y), numel(gT.z));
%
% Nitrogen and silicate, both stored as a field on the [x y z] grid:
%
if any(strcmp(options.fields,'N0'))
    regrid3(gS, gT, pS.pathN0, pT.pathN0, 'N', options.bOverwrite);
end
if any(strcmp(options.fields,'Si0'))
    regrid3(gS, gT, pS.pathSi0, pT.pathSi0, 'Si', options.bOverwrite);
end
%
% The light field, which also has to be interpolated in time:
%
if any(strcmp(options.fields,'parday'))
    if ~isfield(pS,'pathPARday') || ~isfield(pT,'pathPARday')
        fprintf(2,'  parday: %s or %s has no p.pathPARday; skipping\n', ...
            pS.TMname, pT.TMname);
    else
        regridParday(gS, gT, pS, pT, options.bOverwrite);
    end
end
end

% -----------------------------------------------------------------------------
% Interpolate a field on the [x y z] grid
% -----------------------------------------------------------------------------
function regrid3(gS, gT, pathIn, pathOut, sVar, bOverwrite)

if ~check(pathIn, pathOut, sVar, bOverwrite)
    return
end
S = load(pathIn, sVar);
A = S.(sVar);
nNaN = sum(isnan(A(:)));
%
% Fill the gaps over land, level by level, so they are not spread by the
% interpolation:
%
A = fillGaps(gS.x, gS.y, A);

F = griddedInterpolant({gS.x(:), gS.y(:), gS.z(:)}, A, 'linear', 'nearest');
B = F({gT.x(:), gT.y(:), gT.z(:)});
B = max(0, B); % A concentration cannot be negative

report(sVar, A, B, nNaN);
sOut = struct(sVar, B); % so that save stores it under the name sVar
save(pathOut, '-struct', 'sOut');
fprintf('    wrote %s.mat\n', pathOut);
end

% -----------------------------------------------------------------------------
% Interpolate the light field in space and in time
% -----------------------------------------------------------------------------
function regridParday(gS, gT, pS, pT, bOverwrite)

if ~check(pS.pathPARday, pT.pathPARday, 'parday', bOverwrite)
    return
end
S = load(pS.pathPARday, 'parday');
%
% parday is stored as [x y 1 time]; the depth dimension is filled in by
% simulateGlobal. Drop it here and put it back at the end.
%
P = squeeze(S.parday);
nNaN = sum(isnan(P(:)));
P = fillGaps(gS.x, gS.y, P);

nS = size(P,3);
nT = 365/pT.dtTransport;
if abs(nT-round(nT)) > 1e-9
    error('365 days is not a whole number of transport steps for %s.', pT.TMname);
end
nT = round(nT);
%
% The day of the year of each column, with one extra column at either end so
% that the interpolation wraps around new year:
%
tS = (1:nS)'*365/nS;
tS = [tS(1)-365/nS; tS; tS(end)+365/nS];
P  = cat(3, P(:,:,end), P, P(:,:,1));
tT = (1:nT)'*365/nT;

fprintf('    %d time steps of %.4f days -> %d of %.4f days\n', ...
    nS, 365/nS, nT, 365/nT);
F = griddedInterpolant({gS.x(:), gS.y(:), tS}, P, 'linear', 'nearest');
Q = max(0, F({gT.x(:), gT.y(:), tT})); % Light cannot be negative

report('parday', P, Q, nNaN);
% Stored in single precision: it is by far the largest of these files, and the
% light is only used to multiply the uptake rates. simulateGlobal reads it into
% the double array L0, so nothing downstream sees a single.
parday = single(reshape(Q, numel(gT.x), numel(gT.y), 1, nT));
save(pT.pathPARday, 'parday');
fprintf('    wrote %s\n', pT.pathPARday);
end

% -----------------------------------------------------------------------------
% Replace NaNs by the value of the nearest cell that has data. A is [nx ny n],
% where the third dimension is depth or time, and each slice is filled on its
% own. When the pattern of gaps is the same in every slice, which it is for a
% land mask in a time series, the nearest cell is found once and reused.
% -----------------------------------------------------------------------------
function A = fillGaps(x, y, A)

if ~any(isnan(A(:)))
    return
end
[nx, ny, n] = size(A);
[X, Y] = ndgrid(x(:), y(:));
A = reshape(A, nx*ny, n);
mask = isnan(A);

bSame = all(all(mask == mask(:,1)));
if bSame
    iFrom = nearestValid(X, Y, mask(:,1));
    A(mask(:,1),:) = A(iFrom,:);
else
    for k = 1:n
        if ~any(mask(:,k)), continue, end
        iFrom = nearestValid(X, Y, mask(:,k));
        A(mask(:,k),k) = A(iFrom,k);
    end
end
A = reshape(A, nx, ny, n);
end

function iFrom = nearestValid(X, Y, mask)

iValid = find(~mask);
iGap = find(mask);
if isempty(iValid)
    error('A slice of the source field is all NaN.')
end
% Interpolating the indices themselves with 'nearest' gives, for each gap, the
% index of the closest cell that has data:
F = scatteredInterpolant(X(iValid), Y(iValid), iValid, 'nearest', 'nearest');
iFrom = round(F(X(iGap), Y(iGap)));
end

% -----------------------------------------------------------------------------
function bDo = check(pathIn, pathOut, sVar, bOverwrite)

bDo = false;
sIn = pathIn;
if ~endsWith(sIn,'.mat'), sIn = [sIn '.mat']; end
sOut = pathOut;
if ~endsWith(sOut,'.mat'), sOut = [sOut '.mat']; end

if ~isfile(sIn)
    fprintf(2,'  %s: the source file %s is not there; skipping\n', sVar, sIn);
    return
end
if isfile(sOut) && ~bOverwrite
    fprintf('  %s: %s already exists; skipping (use bOverwrite=true)\n', sVar, sOut);
    return
end
fprintf('  %s:\n', sVar);
bDo = true;
end

% -----------------------------------------------------------------------------
function report(sVar, A, B, nNaN)

fprintf('    source: range [%.4g %.4g], mean %.4g, %d NaNs filled (%.1f %%)\n', ...
    min(A(:)), max(A(:)), mean(A(:)), nNaN, 100*nNaN/numel(A));
fprintf('    target: range [%.4g %.4g], mean %.4g, %d NaNs, %.2f %% exact zeros\n', ...
    min(B(:)), max(B(:)), mean(B(:)), sum(isnan(B(:))), 100*mean(B(:)==0));
if sum(isnan(B(:))) > 0
    fprintf(2,'    WARNING: %s still has NaNs\n', sVar);
end
end
