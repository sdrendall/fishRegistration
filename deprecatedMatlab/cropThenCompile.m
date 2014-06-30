function cropThenCompile(startDir, cropCoordPath)

    f = dir(fullfile(startDir, '*.vsi'));
    for i = 1:length(f)
        imPath{i} = fullfile(startDir, f(i).name);
    end
    processedImPaths = cropRemoveZerosAndResize(imPath, cropCoordPath);

    multiTiffBase = fullfile(startDir, 'fullSet_multi');
    populateMultiTiff(processedImPaths, multiTiffBase);