function ai_pt = getAIperTrial( behRes )

Nb = numel( behRes(1).Results );
mvpt = {behRes(1).Results.MaxValuePerTrial};
mvps = cellfun(@max, mvpt) * 1.15;
ai_ptpb = cat(1, mvpt{:} )' ./ mvps;

radAxis = (0:Nb-1)*(2*pi/Nb);
z_axis = exp(1i*radAxis(:));

poly_coords = ai_ptpb .* z_axis.';

ai_pt = arrayfun(@(t) getAIfromCmplxCoords( poly_coords(t,:) ), ...
    1:size( poly_coords, 1 ) );

end

function ai = getAIfromCmplxCoords( z_coords )
ai = area( polyshape( real(z_coords), imag(z_coords) ) );
end