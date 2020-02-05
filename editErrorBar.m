function output = editErrorBar(x, col, width)
    x.mainLine.Color = col;
    x.mainLine.MarkerFaceColor = col;
    x.mainLine.LineWidth = width;
    x.mainLine.MarkerSize = 4;
    x.patch.FaceColor = col;
    x.patch.EdgeColor = 'none';
    x.patch.FaceAlpha = 0.15;
    x.patch.HandleVisibility = 'off';
%     set(x.edge,'Color',x.patch.EdgeColor);
%     x.edge = 'off';
    output = x;
end