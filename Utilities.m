classdef Utilities
    % Utilities  Collection of static utility functions for plotting and colors.
    %
    % Methods:
    %   col = Utilities.getMaterialColor(material)
    %   h   = Utilities.plotByMaterial(x,y,material,Name,Value,...)
    %   h   = Utilities.plotByMaterial(x,y,z,material,Name,Value,...)
    %   Utilities.applyColor(h, material)

    methods (Static)

        function col = getMaterialColor(material)
            if nargin < 1 || isempty(material)
                material = "default";
            end
            m = lower(string(material));
            switch m
                case "cylinder_normal"
                    col = "blue";   
                case "cylinder_rubber"
                    col = "#911eb4";   % purple
                case "cylinder_tpu"
                    col = "#9A6324";
                case "hexagon_normal"
                    col = "blue";
                case "hexagon_rubber"
                    col = "#911eb4";
                case "hexagon_tpu"
                    col = "#9A6324";
                case "oblong_normal"
                    col = "blue";
                case "oblong_rubber"
                    col = "#911eb4"; 
                case "oblong_tpu"
                    col = "#9A6324";
                % case "default"
                %     col = [];                   % let MATLAB choose
                otherwise
                    col = [];                   % unknown -> default
            end
        end

        function h = plotByMaterial(varargin)
            % plotByMaterial  Plot using material color.
            % Usage:
            %   h = Utilities.plotByMaterial(x,y,material,Name,Value,...)
            %   h = Utilities.plotByMaterial(x,y,z,material,Name,Value,...)
            %
            % Detect whether a z vector is provided by argument count and
            % treat the material argument as the last required positional
            % argument before Name-Value pairs.
            narginAll = numel(varargin);
            if narginAll < 3
                error('Not enough input arguments.');
            end

            % Determine if last positional before name-value is numeric vs string
            % We look for first string/char from the end to identify material.
            % Simpler: require that material is provided as a string/char.
            % Find index of material argument (first string/char from left)
            matIdx = [];
            for k = 1:narginAll
                if ischar(varargin{k}) || isstring(varargin{k})
                    matIdx = k;
                    break;
                end
            end
            if isempty(matIdx)
                error('Material name (string) must be supplied as a positional argument.');
            end

            % Extract positional data before material
            posArgs = varargin(1:matIdx-1);
            material = varargin{matIdx};
            nvPairs = {};
            if matIdx < narginAll
                nvPairs = varargin(matIdx+1:end);
            end

            col = Utilities.getMaterialColor(material);

            % Decide 2-D or 3-D based on number of positional args
            if numel(posArgs) == 2
                x = posArgs{1}; y = posArgs{2};
                if isempty(col)
                    h = plot(x, y, nvPairs{:});
                else
                    h = plot(x, y, 'Color', col, nvPairs{:});
                end
            elseif numel(posArgs) == 3
                x = posArgs{1}; y = posArgs{2}; z = posArgs{3};
                if isempty(col)
                    h = plot3(x, y, z, nvPairs{:});
                else
                    h = plot3(x, y, z, 'Color', col, nvPairs{:});
                end
            else
                error('Expected (x,y,material,...) or (x,y,z,material,...)');
            end

            % tag handle(s) with material
            try
                if isgraphics(h)
                    for k = 1:numel(h)
                        h(k).UserData.material = char(string(material));
                    end
                end
            catch
                % ignore if cannot set UserData
            end
        end

        function applyColor(h, material)
            col = Utilities.getMaterialColor(material);
            if isempty(col)
                return
            end

            if isgraphics(h)
                objs = h;
            else
                objs = findall(gcf, '-property', 'Color');
            end

            for k = 1:numel(objs)
                obj = objs(k);
                if isprop(obj, 'FaceColor') && ~isa(obj, 'matlab.graphics.chart.primitive.Line')
                    try
                        obj.FaceColor = col;
                    catch
                        if isprop(obj, 'Color')
                            obj.Color = col;
                        end
                    end
                else
                    if isprop(obj, 'Color')
                        obj.Color = col;
                    end
                end
                try
                    obj.UserData.material = char(string(material));
                catch
                end
            end
        end

    end
end
