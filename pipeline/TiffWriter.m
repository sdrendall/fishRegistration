classdef TiffWriter < handle
    properties
        path
        numIms = 0;
        fileOpen = false;
    end

    methods
        % INIT
        function obj = TiffWriter(tiffPath)
            obj.path = tiffPath;
        end

        % Call to add images
        function append(self, im)
            if self.fileOpen
                self.addToFile(im)
            else
                self.initializeFileWith(im)
            end
        end
    end

    methods (Access=private, Hidden=true)
        function initializeFileWith(self, im)
            imwrite(im, self.path, 'tiff')
            self.fileOpen = true;
            self.numIms = 1;
        end

        function addToFile(self, im)
            imwrite(im, self.path, 'tiff', 'WriteMode', 'append')
            self.numIms = self.numIms + 1;
        end
    end
end