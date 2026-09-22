function [] = correlogram(C, varargin)
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    p = inputParser;
    validX = @(x) validateattributes(x, {'cell', 'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real', '<=', 1, '>=', -1});
    validScalar = @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'scalar','nonempty','real','nonnan'});
    validColorArray = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'nonempty', 'real', 'nonnan', 'ncols', 3});
    addRequired(p, 'C', validX);
    addParameter(p, 'AxisLabels', [], @iscell);
    addParameter(p, 'Labels', [], @iscell);
    addParameter(p, 'cMin', -1, validScalar);
    addParameter(p, 'cMax', 1, validScalar);
    addParameter(p, 'colorMin', [0 0 1], validColorArray);
    addParameter(p, 'colorZero', [1 1 1]*0.96, validColorArray);
    addParameter(p, 'colorMax', [1 0 0], validColorArray);
    addParameter(p, 'numDigits', 2, @isnumeric);
    addParameter(p, 'boxWidthRatio', 1, validScalar);
    addParameter(p, 'boxHeightRatio', 1, validScalar);
    addParameter(p, 'BoxLineWidth', 0.5, validScalar);
    addParameter(p, 'RadiusRatio', 0.9, validScalar);
    addParameter(p, 'Colors', [], @isnumeric);
    addParameter(p, 'VariableNames', cell(length(C), 1), @iscell);
    addParameter(p, 'ShowDiagonal', false, @islogical);
    addParameter(p, 'FontSize', [], validScalar);
    addParameter(p, 'Sorting', true, @islogical);
    addParameter(p, 'XDir', 'normal', @ischar);
    addParameter(p, 'YDir', 'reverse', @ischar);
    parse(p, C, varargin{:});
    param = p.Results;
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));

    [nRow, nColumn] = size(C);
    if((nRow == nColumn) && param.Sorting) 
        for i = 1:2
            Z = linkage(C, 'average');
            D = pdist(C);
            si = optimalleaforder(Z, D);
            C = C(si, si);
            if(~isempty(param.AxisLabels))
                param.AxisLabels = param.AxisLabels(si);
                param.VariableNames = param.VariableNames(si);
            end
        end
    end
    
    if(checkUsingDefaults(p, 'Labels'))
        labels = cell(size(C));
        for iRow = 1:nRow
            for iColumn = 1:nColumn
                cValue = C(iRow, iColumn);
                labels{iRow, iColumn} = num2str(cValue, ['%.', num2str(param.numDigits), 'f']);
                if(iRow == iColumn); labels{iRow, iColumn} = param.VariableNames{iRow}; end
            end
        end
        param.Labels = labels;
    end
   
    if(checkUsingDefaults(p, 'Colors'))
        colorMin = param.colorMin;
        colorZero = param.colorZero;
        colorMax = param.colorMax;
        cols = color_spacing_continuous(C(:), [param.cMin; 0; param.cMax], [colorMin; colorZero; colorMax]);
        param.Colors = reshape(cols, nRow, nColumn, 3);
    end
    
    box_width_ratio = param.boxWidthRatio;
    box_height_ratio = param.boxHeightRatio;
    
    row_height = 1 ./ nRow;
    col_width = 1./ nColumn;
    
    xcPos = zeros(nColumn, 1);
    ycPos = zeros(nRow, 1);
    
    for iRow = 1:nRow
        for iColumn = 1:nColumn
            row_gap = (1 - box_height_ratio) * 0.5;
            col_gap = (1 - box_width_ratio) * 0.5;
            x = 0 + (col_gap + iColumn - 1) * col_width;
            y = 0 + (row_gap + iRow - 1) * row_height;
            w = col_width * (1 - col_gap);
            h = row_height * (1 - row_gap);
            xcenter = x + w/2;
            ycenter = y + h/2;
            xcPos(iColumn) = xcenter;
            ycPos(iRow) = ycenter;
        end
    end

    isHoldOn = ishold();
    if(~isHoldOn)
        cla();
    end

    hold('on');
    set(gca, 'XTick', xcPos, 'XTickLabel', param.AxisLabels);
    set(gca, 'YTick', ycPos, 'YTickLabel', param.AxisLabels);

    xlim([0 1]);
    ylim([0 1]);

    dim1 = get(gcf, 'Position');
    dim2 = get(gca, 'Position');
    axis_width = dim1(3) * dim2(3);
    axis_height = dim1(4) * dim2(4);
    width_to_height = axis_width / axis_height;
    alpha = 0.8;
    if(width_to_height >= 1)
        r_width_mult = 1 / (width_to_height^alpha);
        r_height_mult = 1;
    else
        r_width_mult = 1;
        r_height_mult = (width_to_height^alpha);
    end

    S = struct();
    S.FontSize = 10;
    [txt_width, txt_height] = measureText('-0.99', S, gca);

    scaling = max(txt_width / col_width, txt_height / row_height);
    desired_txt_gap = 0.9;
    font_size = max(floor(S.FontSize * (1/scaling) * desired_txt_gap), 1);

    if(isempty(param.FontSize))
        param.FontSize = font_size;
    end

    txtOptions = struct();
    txtOptions.HorizontalAlignment = 'center';
    txtOptions.VerticalAlignment = 'middle';
    txtOptions.FontSize = param.FontSize;
    
    rec_line_width = param.BoxLineWidth;
    for iRow = 1:nRow
        for iColumn = 1:nColumn
            row_gap = (1 - box_height_ratio) * 0.5;
            col_gap = (1 - box_width_ratio) * 0.5;
            x = 0 + (col_gap + iColumn - 1) * col_width;
            y = 0 + (row_gap + iRow - 1) * row_height;
            w = col_width * (1 - col_gap);
            h = row_height * (1 - row_gap);
            xcenter = x + w/2;
            ycenter = y + h/2;
            xcPos(iColumn) = xcenter;
            ycPos(iRow) = ycenter;
            pos = [x y w h];
            color = param.Colors(iRow, iColumn, :);
            
            if(iRow > iColumn)
                rectangle('Position',pos, 'FaceColor', color, 'LineWidth', rec_line_width)
                txt = param.Labels{iRow, iColumn};
                text(xcenter, ycenter, txt, txtOptions);
            end
            if((iRow == iColumn) && param.ShowDiagonal)
                rectangle('Position',pos, 'FaceColor', color, 'LineWidth', rec_line_width)
                txt = param.Labels{iRow, iColumn};
                text(xcenter, ycenter, txt, txtOptions);
            end
            if(iRow < iColumn)
                cvalue = C(iRow, iColumn);
                theta_start = -0.5 * pi;
                theta_end = theta_start + (cvalue * 2 * pi);
                theta_range = theta_end - theta_start;
                
                n = 100;
                ww = w * param.RadiusRatio * r_width_mult;
                hh = h * param.RadiusRatio * r_height_mult;
                
                theta = [(0:n-1)]*(2*pi)/n;
                xx = xcenter + [0.5*ww*cos(theta)];
                yy = ycenter + [0.5*hh*sin(theta)];
                Pcircle = polyshape(xx,yy);

                theta = [(0:n-1) (n-1e-3)]*(theta_range)/n + theta_start;
                xx = xcenter + [0 0.5*ww*cos(theta)];
                yy = ycenter + [0 0.5*hh*sin(theta)];
                P = polyshape(xx,yy);
                
                plot(Pcircle, 'FaceAlpha', 0, 'LineStyle', '-', 'EdgeColor', [1 1 1] * 0.55);
                plot(P, 'FaceAlpha', 1, 'FaceColor', color);
            end
            
        end
    end
    set(gca, 'XDir', param.XDir);
    set(gca, 'YDir', param.YDir)
    if(~isHoldOn); hold('off'); end
end

