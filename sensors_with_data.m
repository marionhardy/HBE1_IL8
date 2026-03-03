function keep = sensors_with_data(Tcond, cellLine, sensorVars)
Tc = Tcond(Tcond.cell_line == cellLine, :);
keep = sensorVars;
keep = keep(arrayfun(@(v) any(~isnan(Tc.(v))), keep));
end