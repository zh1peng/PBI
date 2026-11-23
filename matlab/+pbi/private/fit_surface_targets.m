function fit_surface_targets(SPM, PBI, outputDir, contrasts, bold, meshes)
% Use the same SPM design for each target geometry.
sides = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
targetBase = {pwd};
isCifti = isscalar(bold) && endsWith(bold, '.dtseries.nii');
if isCifti
    assert(exist('cifti_read', 'file') == 2, 'Add Washington University cifti-matlab to the MATLAB path.');
    dt = cifti_read(char(bold));
    dense = dt.diminfo{1};
    time = dt.diminfo{2};
    assert(strcmp(dense.type,'dense') && strcmp(time.type,'series') && ...
        strcmp(time.seriesUnit,'SECOND') && time.length == SPM.nscan && ...
        abs(time.seriesStep - SPM.xY.RT) < 1e-5, 'CIFTI must match the source scan count and TR.');
    [surfaces, volumes] = cifti_diminfo_dense_get_structures(dense);
    assert(all(ismember(surfaces, sides)), 'Supported CIFTI surfaces are CORTEX_LEFT and CORTEX_RIGHT.');
    maps = nan(dense.length, numel(contrasts), 2, 'single');
else
    assert(numel(bold) == 2 && all(endsWith(bold(:), '.gii')), ...
        'Supply two GIFTI time series (left, right) or one CIFTI .dtseries.nii.');
    surfaces = sides;
end

for h = 1:numel(sides)
    side = sides{h};
    if ~ismember(side, surfaces), continue; end
    assert(numel(meshes) == 2, 'Supply TargetMeshes in [left, right] order.');
    meshFile = spm_select('CPath', char(meshes(h)), targetBase);
    mesh = gifti(meshFile);
    nVertices = size(mesh.vertices, 1);
    if isCifti
        info = cifti_diminfo_dense_get_surface_info(dense, side);
        assert(info.numverts == nVertices, 'Mesh vertex count differs from CIFTI for %s.', side);
        vertices = info.vertlist1;
        values = dt.cdata(info.ciftilist, :);
    else
        series = gifti(spm_select('CPath', char(bold(h)), targetBase));
        values = series.cdata;
        if isfield(series, 'indices')
            vertices = series.indices;
        else
            vertices = 1:nVertices;
        end
    end
    assert(isequal(size(values), [numel(vertices), SPM.nscan]) && ...
        all(vertices >= 1 & vertices <= nVertices) && numel(unique(vertices)) == numel(vertices), ...
        'Vertex mapping or scan count differs for %s.', side);
    folder = fullfile(outputDir, side);
    scans = surface_frames(folder, meshFile, nVertices, vertices, values);
    PBI_1st_level(SPM, PBI, folder, contrasts, scans);
    if isCifti && ~isempty(contrasts)
        maps(info.ciftilist,:,:) = read_maps(folder, '.gii', vertices, numel(contrasts));
    end
end

if ~isCifti, return; end
if ~isempty(volumes)
    folder = fullfile(outputDir, 'volume');
    info = cifti_diminfo_dense_get_volume_all_info(dense, true);
    voxels = cifti_vox2ind(info.voldims, info.voxlist1);
    scans = volume_frames(folder, info, voxels, dt.cdata);
    PBI_1st_level(SPM, PBI, folder, contrasts, scans);
    if ~isempty(contrasts)
        maps(info.ciftilist,:,:) = read_maps(folder, '.nii', voxels, numel(contrasts));
    end
end
prefix = {'con', 'spmT'};
for c = 1:numel(contrasts)
    for k = 1:numel(prefix)
        file = fullfile(outputDir, sprintf('%s_%04d.dscalar.nii', prefix{k}, c));
        cifti_write_from_template(dt, maps(:,c,k), file, 'namelist', {contrasts(c).name});
    end
end
end

function scans = surface_frames(folder, meshFile, nVertices, vertices, values)
frameDir = fullfile(folder, 'frames');
mkdir(frameDir);
frame = nan(nVertices, 1, 'single');
scans = cell(size(values,2), 1);
for t = 1:numel(scans)
    frame(vertices) = values(:,t);
    g = gifti(frame);
    g.private.metadata(1).name = 'SurfaceID';
    g.private.metadata(1).value = meshFile;
    scans{t} = fullfile(frameDir, sprintf('frame_%04d.gii', t));
    save(g, scans{t}, 'ExternalFileBinary');
end
end

function scans = volume_frames(folder, info, voxels, data)
mkdir(folder);
header = struct('fname', fullfile(folder, 'input.nii'), 'dim', info.voldims, ...
    'dt', [spm_type('float32') 0], 'pinfo', [1;0;0], 'mat', info.volsform1);
frame = zeros(info.voldims, 'single');
scans = cell(size(data,2), 1);
for t = 1:numel(scans)
    frame(voxels) = data(info.ciftilist,t);
    header.n = [t 1];
    spm_write_vol(header, frame);
    scans{t} = sprintf('%s,%d', header.fname, t);
end
end

function maps = read_maps(folder, ext, indices, nContrast)
% Leave unestimated grayordinates missing in the merged maps.
mask = spm_data_read(spm_data_hdr_read(fullfile(folder, ['mask' ext])), indices) > 0;
maps = nan(numel(indices), nContrast, 2, 'single');
prefix = {'con', 'spmT'};
for c = 1:nContrast
    for k = 1:numel(prefix)
        file = fullfile(folder, sprintf('%s_%04d%s', prefix{k}, c, ext));
        values = spm_data_read(spm_data_hdr_read(file), indices);
        maps(mask,c,k) = values(mask);
    end
end
end
