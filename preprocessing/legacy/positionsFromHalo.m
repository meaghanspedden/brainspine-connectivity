function [positions, goodPositions, distError, distances, angleError, angles] = positionsFromHalo(cfg)
% Reads in and (optionally) aligns halo position information before
% exporting as a positions.tsv file compatible with spm_opm_create.
% FORMAT [positions, goodPositions, distError, distances, angleError, angles] = positionsFromHalo(cfg)
%
% Input Parameters:
%     cfg:				Struct
%         cfg.haloResults:	string, full path to RESULTS_Localization Solver*.csv.
%         cfg.output:		string, full path to output for positions.tsv.
%         cfg.align:		boolean, attempt to align to anatomy? default = false.
%		  cfg.alignMethod:	string, 'Kabsch' [default], 'RotationOptimised',
%							'TranslationOptimised', 'TransformOptimised'.
%         cfg.cadPositions:	string, full path to existing positions.tsv from CAD
%         cfg.plot:			boolean, plot channel frames. default = false
% 							If cfg.align then outliers are shown in red and 
%							good positions in blue. If cadPositions are provided
% 							they are shown in black. 
%		  cfg.threshold		int, threshold to define outliers. default = 10
%		  cfg.mod			int, inverse of how many pairs have to exceed threshold. 
%							default = 6.
%		  cfg.lvm			string, full path to .lvm meg recording. This is necessary
%							to align channel names with different naming conventions 
%							used by Quspin. If it is not provided then a guess at channel names
%							is made.
%		  cfg.checkFlip		boolean, default = false. Checks if halo orientations
%							are flipped
% %TODO	  cfg.anatomy		string, full path to anatomical information. 
% Output:
%     positions: A table following the format used in positions.tsv
%     compatible with spm_opm_create. 
%     goodPositions: Empty unless cfg.align = true. Otherwise it contains
%     a subset of positions, with outliers removed. Outliers may be due to
%     the halo software being unable to localise the sensor successfully,
%     or producing results vastly different to cad. If you are getting
%     strange results chance cfg.mod/cfg.theshold. 
%     distError: mean distance error between HALO and CAD.
%	  distances: distances indexed to goodPositions rows.
%	  angleError: same as distError, but angles in degrees.
%	  angles: same as distances, but angles in degrees.
%__________________________________________________________________________
%
% Further help:
% positionsFromHalo is designed to integrate the calibration table from the
% HALO system (https://doi.org/10.48550/arXiv.2410.08718) with existing
% pipelines that had used CAD, previously aligned with anatomy, for spatial
% coregistration. In particular, if you have both CAD information (labelled
% or unlabelled) this function can be used to get the correct
% position/orientation information provided by HALO while maintaining
% coregistation to anatomy given by CAD. Note, the alignment is
% approximate.
%_________________________________________________________________________

% Nicholas Alexander (n.alexander@ucl.ac.uk)
% Adapted code from info2pos_neuro1 written by Stephanie Mellor
% Copyright (C) 2024 Department of Imaging Neuroscience, UCL

%==========================================================================
% - P O S I T I O N S   F R O M   H A L O
%==========================================================================

%-CHECK function inputs
%--------------------------------------------------------------------------
if ~isfield(cfg,'plot')
	cfg.plot = false;
end
if ~isfield(cfg,'align')
	cfg.align = false;
end
if ~isfield(cfg,'threshold')
	cfg.threshold = 10;
end
if isfield(cfg,'cadPositions') && ~isempty(cfg.cadPositions)
	cfg.align = true;
end
if ~isfield(cfg,'output') || isempty(cfg.output)
	warning('Setting output to positions.tsv in cd. Please specify cfg.output.')
	cfg.output = fullfile(cd,'positions.tsv');
elseif ~strcmp(cfg.output(end-3:end),'.tsv')
	disp('Setting output to end in .tsv')
	cfg.output = [cfg.output, '.tsv'];
end
if ~isfield(cfg,'haloResults') || isempty(cfg.haloResults)
	error(['Please provide cfg.haloResults. e.g. ' ...
		'..\Results 2024-10-10 11-20\RESULTS_Localization Solver-2024-10-10 11-20.csv'])
end
if ~(isfield(cfg,'cadPositions') && ~isempty(cfg.cadPositions)) && ...
		~(isfield(cfg,'anatomy') && ~isempty(cfg.anatomy))
	cfg.align = false;
else
	if isa(cfg.cadPositions,'string') || isa(cfg.cadPositions,'char')
		cadPositions = readtable(cfg.cadPositions, "FileType","text",'Delimiter', '\t');
		cfg.cadPositions = cadPositions;
	else
		cadPositions = cfg.cadPositions;
	end
end
if isfield(cfg,'align') && (~isfield(cfg,'alignMethod') || isempty(cfg.alignMethod))
	cfg.alignMethod = 'Kabsch';
	disp('Defaulting to Kabsch alignment')
end
if isfield(cfg,'alignMethod') && strcmp(cfg.alignMethod, 'Optimised')
	warning('Optimised is no longer an option. Defaulting to TransformOptimised');
	cfg.alignMethod = 'TransformOptimised';
end
if ~isfield(cfg,'mod') || isempty(cfg.mod)
	cfg.mod = 6;
end
if ~isfield(cfg,'lvm')
	cfg.lvm = [];
end
if ~isfield(cfg,'checkFlip') || isempty(cfg.checkFlip)
	cfg.checkFlip = true;
end



%-PROCESS halo results file
%--------------------------------------------------------------------------
haloResults = readtable(cfg.haloResults, "FileType","text",'Delimiter', ...
				'comma','VariableNamingRule','preserve');
axesLabels = {'X', 'Y', 'Z'};
numAxes = length(axesLabels);
numObjects = length(haloResults.Object);

positions = struct('name', cell(1, numObjects * numAxes), ...
                   'Px', zeros(1, numObjects * numAxes), ...
                   'Py', zeros(1, numObjects * numAxes), ...
                   'Pz', zeros(1, numObjects * numAxes), ...
                   'Ox', zeros(1, numObjects * numAxes), ...
                   'Oy', zeros(1, numObjects * numAxes), ...
                   'Oz', zeros(1, numObjects * numAxes));
badHaloIdx = zeros(1, numObjects * numAxes);

for objIdx = 1:numObjects
    % Convert positions to correct units
    posX = haloResults.x_pos(objIdx) * 1000;
    posY = haloResults.y_pos(objIdx) * 1000;
    posZ = haloResults.z_pos(objIdx) * 1000;

    % Assign orientation matrix values
    orientationMatrix = [
        haloResults.Nxx(objIdx), haloResults.Nxy(objIdx), haloResults.Nxz(objIdx);
        haloResults.Nyx(objIdx), haloResults.Nyy(objIdx), haloResults.Nyz(objIdx);
        haloResults.Nzx(objIdx), haloResults.Nzy(objIdx), haloResults.Nzz(objIdx)
    ];

    for axisIdx = 1:numAxes
        % Determine current struct index
        newObjIdx = (objIdx - 1) * numAxes + axisIdx;

        % Set combined name and positions
        positions(newObjIdx).name = [haloResults.Object{objIdx}, '-', axesLabels{axisIdx}];
        positions(newObjIdx).Px = posX;
        positions(newObjIdx).Py = posY;
        positions(newObjIdx).Pz = posZ;

        % Extract and normalise orientation vector for current axis
        orientationVector = orientationMatrix(axisIdx, :);
        normFactor = norm(orientationVector);
        positions(newObjIdx).Ox = orientationVector(1) / normFactor;
        positions(newObjIdx).Oy = orientationVector(2) / normFactor;
        positions(newObjIdx).Oz = orientationVector(3) / normFactor;

        % Assign halo quality score
        badHaloIdx(newObjIdx) = haloResults.obj_func_value(objIdx);
    end
end

% 0 seems to code for failed localisation
badHaloIdx = (badHaloIdx == 0);
positions = struct2table(positions);

%-ALIGN if requested
%--------------------------------------------------------------------------
if (cfg.align)
	goodPositions = positions(~badHaloIdx,:);

	if isfield(cfg,'cadPositions') && ~isempty(cfg.cadPositions)
		% Read in and format the cadPositions		
		cadPositions.name = cellfun(@(x) x(end-3:end), cadPositions.name, 'UniformOutput', false);

		% Match index
		[~,b] = ismember(goodPositions.name,cadPositions.name);

		if ~isempty(b)
			cadPositions = cadPositions(b(b ~= 0),:);
			goodPositions = goodPositions(b ~= 0, :);
		else
			% Get relative distances between channels
			cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
			haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
			cadDists = pdist(cadPos);
			haloDists = pdist(haloPos);
			
			% Sort the distances and find closest matches
			[~, idxCad] = sort(cadDists);
			[~, idxHalo] = sort(haloDists);
			maxMatches = min(length(idxCad), length(idxHalo));  % Number of pairs to match
			matchedIdxCad = idxCad(1:maxMatches);  % Indices of closest pairs in A
			matchedIdxHalo = idxHalo(1:maxMatches);  % Indices of closest pairs in B
			
			% Identify pairs
			[rowCad, colCad] = ind2sub(size(cadPos, 1), matchedIdxCad);
			[rowHalo, colHalo] = ind2sub(size(haloPos, 1), matchedIdxHalo);
			pointsCad = unique([rowCad; colCad]);
			pointsHalo = unique([rowHalo; colHalo]);
			
			% Update table
			numMatches = min(length(pointsCad), length(pointsHalo));
			cadPositions = cadPositions(pointsCad(1:numMatches), :);
			goodPositions = goodPositions(pointsHalo(1:numMatches), :);
		end
		
		% Get relative distances between channels
		cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
		haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
		cadDists = pdist(cadPos);
		haloDists = pdist(haloPos);
		diffDists = cadDists - haloDists;

		% Identify shared and outlier positions
		sharedPosIdx = (abs(diffDists) == 0);
		if ~all(sharedPosIdx)
			outlierIdx = (abs(diffDists) > cfg.threshold);
			sharedPairs = squareform(sharedPosIdx);
			outlierPairs = squareform(outlierIdx);
			outlierPairs = ~sharedPairs + ~outlierPairs;
			maxPairs = (length(outlierPairs) + 1) * 2;
			relativeOutlierCount = -sum(outlierPairs,1) + maxPairs;
			newBadHaloIdx = relativeOutlierCount > maxPairs / cfg.mod;
			goodPositions = goodPositions(~newBadHaloIdx,:);
		else
			warning('It is very likely you are comparing two of the same positions files.')
		end
		

		% Match index (again)
		[~,b] = ismember(goodPositions.name,cadPositions.name);

		% Depends on whether slot2sens has been recorded and applied
		if ~isempty(b)
			cadPositions = cadPositions(b,:);
		else
			% Get relative distances between channels
			cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
			haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
			cadDists = pdist(cadPos);
			haloDists = pdist(haloPos);
			
			% Sort the distances and find closest matches
			[~, idxCad] = sort(cadDists);
			[~, idxHalo] = sort(haloDists);
			maxMatches = min(length(idxCad), length(idxHalo));  % Number of pairs to match
			matchedIdxCad = idxCad(1:maxMatches);  % Indices of closest pairs in A
			matchedIdxHalo = idxHalo(1:maxMatches);  % Indices of closest pairs in B
			
			% Identify pairs
			[rowCad, colCad] = ind2sub(size(cadPos, 1), matchedIdxCad);
			[rowHalo, colHalo] = ind2sub(size(haloPos, 1), matchedIdxHalo);
			pointsCad = unique([rowCad; colCad]);
			pointsHalo = unique([rowHalo; colHalo]);
			
			% Update table
			numMatches = min(length(pointsCad), length(pointsHalo));
			cadPositions = cadPositions(pointsCad(1:numMatches), :);
			goodPositions = goodPositions(pointsHalo(1:numMatches), :);
		end

		% Subtract centroid
		haloCent(1) = mean(goodPositions.Px);
		haloCent(2) = mean(goodPositions.Py);
		haloCent(3) = mean(goodPositions.Pz);
		goodPositions.Px = goodPositions.Px - haloCent(1);
		goodPositions.Py = goodPositions.Py - haloCent(2);
		goodPositions.Pz = goodPositions.Pz - haloCent(3);
		positions.Px = positions.Px - haloCent(1);
		positions.Py = positions.Py - haloCent(2);
		positions.Pz = positions.Pz - haloCent(3);
		
		cadCent(1) = mean(cadPositions.Px);
		cadCent(2) = mean(cadPositions.Py);
		cadCent(3) = mean(cadPositions.Pz);
		cadPositions.Px = cadPositions.Px - cadCent(1);
		cadPositions.Py = cadPositions.Py - cadCent(2);
		cadPositions.Pz = cadPositions.Pz - cadCent(3);
		
		% Align
		switch cfg.alignMethod
			case 'Kabsch'
				cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
				haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
				Boriginal = [positions.Px, positions.Py, positions.Pz];

				C = cadPos' * haloPos;
				[U, ~, V] = svd(C);
				R = V * U';
				if det(R) < 0
					V(:,3) = -V(:,3);
					R = V * U';
				end

			case 'RotationOptimised' % Previously 'Optimised'
    			% Initial alignment
    			cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
    			haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
    			Boriginal = [positions.Px, positions.Py, positions.Pz];

    			C = cadPos' * haloPos;
    			[U, ~, V] = svd(C);
    			R = V * U';
    			if det(R) < 0
        			V(:,3) = -V(:,3);
        			R = V * U';
    			end
    			
				% Optimise rotation to minimise error
    			initialRotationParams = [0, 0, 0];
    			objectiveFunc = @(rotationParams) errorFromRotate(rotationParams, haloPos, cadPos, R);
    			options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
    			optimalRotationParams = fminunc(objectiveFunc, initialRotationParams, options);
			
    			% Apply optimal rotation
    			thetaX = optimalRotationParams(1);
    			thetaY = optimalRotationParams(2);
    			thetaZ = optimalRotationParams(3);
			
    			Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
    			Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
    			Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];
    			R = Rx * Ry * Rz * R;

			case 'TranslationOptimised'
    			% Initial alignment
    			cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
    			haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
    			Boriginal = [positions.Px, positions.Py, positions.Pz];

    			C = cadPos' * haloPos;
    			[U, ~, V] = svd(C);
    			R = V * U';
    			if det(R) < 0
        			V(:,3) = -V(:,3);
        			R = V * U';
    			end
    			
				haloPos = haloPos * R;
				Boriginal = Boriginal * R;

				% Optimise translation to minimise error
    			initialTranslationParams = [0, 0, 0];
    			objectiveFunc = @(translationParams) errorFromTranslate(translationParams, haloPos, cadPos);
    			options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
    			optimalTranslationParams = fminunc(objectiveFunc, initialTranslationParams, options);
			
    			 % Translate
    			xTran = optimalTranslationParams(1);
    			yTran = optimalTranslationParams(2);
    			zTran = optimalTranslationParams(3);
			
				haloPos(:,1) = haloPos(:,1) + xTran;
				haloPos(:,2) = haloPos(:,2) + yTran;
				haloPos(:,3) = haloPos(:,3) + zTran;

				Boriginal(:,1) = Boriginal(:,1) + xTran;
				Boriginal(:,2) = Boriginal(:,2) + yTran;
				Boriginal(:,3) = Boriginal(:,3) + zTran;

				% Final alignment
				C = cadPos' * haloPos;
				[U, ~, V] = svd(C);
				R = V * U';
				if det(R) < 0
					V(:,3) = -V(:,3);
					R = V * U';
				end
			
			case 'TransformOptimised'
    			% Initial alignment
    			cadPos = [cadPositions.Px, cadPositions.Py, cadPositions.Pz];
    			haloPos = [goodPositions.Px, goodPositions.Py, goodPositions.Pz];
    			Boriginal = [positions.Px, positions.Py, positions.Pz];

    			C = cadPos' * haloPos;
    			[U, ~, V] = svd(C);
    			R = V * U';
    			if det(R) < 0
        			V(:,3) = -V(:,3);
        			R = V * U';
    			end
    			
				haloPos = haloPos * R;
				Boriginal = Boriginal * R;

				% Optimise transform to minimise error
    			initialTransformParams = [0, 0, 0, 0, 0, 0];
    			objectiveFunc = @(translationParams) errorFromTransform(translationParams, haloPos, cadPos);
    			options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
    			optimalTransformParams = fminunc(objectiveFunc, initialTransformParams, options);
			
    			% Translate
    			xTran = optimalTransformParams(1);
    			yTran = optimalTransformParams(2);
    			zTran = optimalTransformParams(3);
			
				haloPos(:,1) = haloPos(:,1) + xTran;
				haloPos(:,2) = haloPos(:,2) + yTran;
				haloPos(:,3) = haloPos(:,3) + zTran;

				Boriginal(:,1) = Boriginal(:,1) + xTran;
				Boriginal(:,2) = Boriginal(:,2) + yTran;
				Boriginal(:,3) = Boriginal(:,3) + zTran;

				% Rotate
    			thetaX = optimalTransformParams(4);
    			thetaY = optimalTransformParams(5);
    			thetaZ = optimalTransformParams(6);
    			Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
    			Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
    			Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];
			
				% Second alignment
				C = cadPos' * haloPos;
				[U, ~, V] = svd(C);
				R = V * U';
				if det(R) < 0
					V(:,3) = -V(:,3);
					R = V * U';
				end
			
    			% Rotate haloPos and compute error
    			R = Rx * Ry * Rz * R;

			otherwise
				error('Check cfg.alignMethod');

		end
		
		haloPos = haloPos * R;
		Boriginal = Boriginal * R;
		
		% Report error
		distances = sqrt(sum((haloPos - cadPos).^2, 2));
		distError = mean(distances);
		disp(['Mean error between CAD and Halo positions: ', num2str(distError)]);

		% Update halo positions
		goodPositions.Px = haloPos(:,1) + cadCent(1);
		goodPositions.Py = haloPos(:,2) + cadCent(2);
		goodPositions.Pz = haloPos(:,3) + cadCent(3);
		C = [goodPositions.Ox, goodPositions.Oy, goodPositions.Oz];
		C = C * R;
		goodPositions.Ox = C(:,1);
		goodPositions.Oy = C(:,2);
		goodPositions.Oz = C(:,3);

		% And update the original as well
		positions.Px = Boriginal(:,1) + cadCent(1);
		positions.Py = Boriginal(:,2) + cadCent(2);
		positions.Pz = Boriginal(:,3) + cadCent(3);
		C = [positions.Ox, positions.Oy, positions.Oz];
		C = C * R;
		positions.Ox = C(:,1);
		positions.Oy = C(:,2);
		positions.Oz = C(:,3);
	
		% Reapply offset to cad
		cadPositions.Px = cadPos(:,1) + cadCent(1);
		cadPositions.Py = cadPos(:,2) + cadCent(2);
		cadPositions.Pz = cadPos(:,3) + cadCent(3);

		% Measure orientation error
		U = [cadPositions.Ox, cadPositions.Oy, cadPositions.Oz];
		V = [goodPositions.Ox, goodPositions.Oy, goodPositions.Oz];
		
		dotProducts = sum(U .* V, 2); 
		normU = sqrt(sum(U.^2, 2));
		normV = sqrt(sum(V.^2, 2));
		
		% Compute cosine of the angle between each pair of vectors
		cosTheta = dotProducts ./ (normU .* normV);
		cosTheta = max(min(cosTheta, 1), -1);
		
		% Compute the angles (in radians) between the vectors
		angleRadians = acos(cosTheta);
		angles = rad2deg(angleRadians);
		idxX = endsWith(goodPositions.name, '-X');
		idxY = endsWith(goodPositions.name, '-Y');
		idxZ = endsWith(goodPositions.name, '-Z');

		angleError(1) = mean(angles(idxX));
		angleError(2) = mean(angles(idxY));
		angleError(3) = mean(angles(idxZ));
		disp(['Mean error between CAD and Halo orientations: x, ', num2str(angleError(1)),...
				'; y, ', num2str(angleError(2)), '; z, ', num2str(angleError(3))]);
		
		if cfg.checkFlip
    		flipOptions = [1, 1, 1;...
							-1, 1, 1;...
							1, -1, 1;...
							1, 1, -1;...
							-1, -1, 1;...
							-1, 1, -1;...
							1, -1, -1;...
							-1, -1, -1];
			flipDesc = {'No flip';'Flip Vx';'Flip Vy';'Flip Vz';
						'Flip Vx and Vy';'Flip Vx and Vz';'Flip Vy and Vz';
						'Flip all axes'};
		
    		bestAngleError = angleError; % Initialise with the original error
    		bestFlip = [1, 1, 1]; % Track the best flip configuration
    		bestFlipDesc = 1;

    		% Iterate over all flip combinations
    		for f = 1:size(flipOptions, 1)
        		flip = flipOptions(f, :); % Current axis flip configuration
	
        		% Apply permutation and flip
        		Vflipped = V .* flip;
	
        		% Compute dot products and vector norms
        		dotProducts = sum(U .* Vflipped, 2);
        		normU = sqrt(sum(U.^2, 2));
        		normV = sqrt(sum(Vflipped.^2, 2));
	
        		% Compute cosine of the angle between each pair of vectors
        		cosTheta = dotProducts ./ (normU .* normV);
        		cosTheta = max(min(cosTheta, 1), -1);
	
        		% Compute the angles (in radians) and their error
        		angleRadians = acos(cosTheta);
        		anglesFlipped = rad2deg(angleRadians);
        		angleErrorFlipped = mean(anglesFlipped);
	
        		% Check if this configuration yields a smaller error
        		if angleErrorFlipped < bestAngleError
            		bestAngleError = angleErrorFlipped;
            		bestFlip = flip;
					bestFlipDesc = f;
        		end
    		end
		
    		disp(['Mean error between CAD and Halo is lowest with flip condition: ',flipDesc{bestFlipDesc},'; at ', num2str(bestAngleError)]);
		
    		% Apply the best configuration if it improves the error
    		flipped = any(bestFlip ~= [1, 1, 1]);
    		if flipped
        		warning('Halo positions are likely flipped, or CAD is incorrect.');
        		warning('Returned output has been modified.');
		
        		goodPositions.Ox = -goodPositions.Ox;
				goodPositions.Oy = -goodPositions.Oy;
				goodPositions.Oz = -goodPositions.Oz;
				positions.Ox = -positions.Ox;
				positions.Oy = -positions.Oy;
				positions.Oz = -positions.Oz;

    		end
		else
    		flipped = false;
		end

		% Read in original CAD
		originalCadPositions = cfg.cadPositions;

		% Visualise
		if cfg.plot
			if flipped
				warning('Plotting flipped halo positions')
			end
			figure
			hold on
			scale = 10; % TODO: scale with unit
			quiver3(originalCadPositions.Px, originalCadPositions.Py, originalCadPositions.Pz,...
				originalCadPositions.Ox*scale, originalCadPositions.Oy*scale , originalCadPositions.Oz*scale,...
				'off','Color','g');
			quiver3(cadPositions.Px, cadPositions.Py, cadPositions.Pz,...
				cadPositions.Ox*scale, cadPositions.Oy*scale , cadPositions.Oz*scale,...
				'off','Color','k');
			quiver3(positions.Px, positions.Py, positions.Pz,...
				positions.Ox*scale, positions.Oy*scale, positions.Oz*scale,...
				'off','Color','r');
			quiver3(goodPositions.Px, goodPositions.Py, goodPositions.Pz,...
				goodPositions.Ox*scale, goodPositions.Oy*scale, goodPositions.Oz*scale,...
				'off','Color','b');
			legend({'Bad CAD', 'Good CAD','Bad Halo','Good Halo'})
			axis equal
			hold off
		end

		% Add channel number info from lvm file.
		if ~isempty(cfg.lvm)
			positions = getchannelNumberFromLvm(cfg.lvm,positions);
			goodPositions = getchannelNumberFromLvm(cfg.lvm,goodPositions);
		else
			warning('Guessing channel number based on A-H, 1-8 format. Supply cfg.lvm to avoid this.')
			for posIdx = 1:length(positions.name)
				boardNameNumber = positions(posIdx,:).name{1}(1:2);
				positionNumber = getChannelNumberFromLvm(boardNameNumber);
				positions(posIdx,:).name{1} = [num2str(positionNumber),'-',positions(posIdx,:).name{1}];
			end
			for posIdx = 1:length(goodPositions.name)
				boardNameNumber = goodPositions(posIdx,:).name{1}(1:2);
				positionNumber = getChannelNumberFromLvm(boardNameNumber);
				goodPositions(posIdx,:).name{1} = [num2str(positionNumber),'-',goodPositions(posIdx,:).name{1}];
			end

		end
	elseif isfield(cfg,'anatomy') && ~isempty(cfg.anatomy)
		% To do...
	end
else
	if ~isempty(cfg.lvm)
		positions = getchannelNumberFromLvm(lvm,positions);
	else
		warning('Guessing channel number based on A-H, 1-8 format. Supply cfg.lvm to avoid this.')
		for posIdx = 1:length(positions.name)
			boardNameNumber = positions(posIdx,:).name{1}(1:2);
			positionNumber = getChannelNumberFromLvm(boardNameNumber);
			positions(posIdx,:).name{1} = [num2str(positionNumber),'-',positions(posIdx,:).name{1}];
		end
	end
	% Empty outputs if no alignment
	goodPositions = [];
	distError = [];
	distances = [];
	angles = [];
	angleError  = [];
	% Visualise
	if cfg.plot
		figure
		hold on
		scale = 10;
		quiver3(positions.Px, positions.Py, positions.Pz,...
			positions.Ox*scale, positions.Oy*scale, positions.Oz*scale,...
			'off','Color','r');
		
		axis equal
		hold off
	end
end


%-SAVE output/s
%--------------------------------------------------------------------------
writetable(positions,cfg.output,'Delimiter','tab','QuoteStrings',...
    false,'FileType', 'text');
if cfg.align
	writetable(goodPositions,[cfg.output(1:end-4),'_good.tsv'],'Delimiter','tab','QuoteStrings',...
    false,'FileType', 'text');
end
end

%==========================================================================
% - E R R O R   F R O M   R O T A T E
%==========================================================================
function error = errorFromRotate(rotationParams, haloPos, cadPos, R)
    % Rotate
    thetaX = rotationParams(1);
    thetaY = rotationParams(2);
    thetaZ = rotationParams(3);
    Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
    Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
    Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];
    Ropt = Rx * Ry * Rz * R;
    transformedHaloPos = haloPos * Ropt;

	% Compute error
    distances = sqrt(sum((transformedHaloPos - cadPos).^2, 2));
    error = mean(distances);
end

%==========================================================================
% - E R R O R   F R O M   T R A N S L A T E
%==========================================================================
function error = errorFromTranslate(translationParams, haloPos, cadPos)
    % Translate
    xTran = translationParams(1);
    yTran = translationParams(2);
    zTran = translationParams(3);

	haloPos(:,1) = haloPos(:,1) + xTran;
	haloPos(:,2) = haloPos(:,2) + yTran;
	haloPos(:,3) = haloPos(:,3) + zTran;

	% Second alignment
	C = cadPos' * haloPos;
	[U, ~, V] = svd(C);
	R = V * U';
	if det(R) < 0
		V(:,3) = -V(:,3);
		R = V * U';
	end
    transformedHaloPos = haloPos * R;

	% Compute error
    distances = sqrt(sum((transformedHaloPos - cadPos).^2, 2));
    error = mean(distances);
end

%==========================================================================
% - E R R O R   F R O M   T R A N S F O R M
%==========================================================================
function error = errorFromTransform(transformParams, haloPos, cadPos)
    % Translate
    xTran = transformParams(1);
    yTran = transformParams(2);
    zTran = transformParams(3);

	haloPos(:,1) = haloPos(:,1) + xTran;
	haloPos(:,2) = haloPos(:,2) + yTran;
	haloPos(:,3) = haloPos(:,3) + zTran;

	% Rotate
    thetaX = transformParams(4);
    thetaY = transformParams(5);
    thetaZ = transformParams(6);
    Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
    Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
    Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];

	% Second alignment
	C = cadPos' * haloPos;
	[U, ~, V] = svd(C);
	R = V * U';
	if det(R) < 0
		V(:,3) = -V(:,3);
		R = V * U';
	end
    Ropt = Rx * Ry * Rz * R;
    transformedHaloPos = haloPos * Ropt;

	% Compute error
    distances = sqrt(sum((transformedHaloPos - cadPos).^2, 2));
    error = mean(distances);
end
% 
% %==========================================================================
% % - E R R O R   F R O M   T R A N S F O R M
% %==========================================================================
% function error = errorFromTransformScale(transformParams, haloPos, cadPos)
%     % Translate
%     xTran = transformParams(1);
%     yTran = transformParams(2);
%     zTran = transformParams(3);
% 
% 	haloPos(:,1) = haloPos(:,1) + xTran;
% 	haloPos(:,2) = haloPos(:,2) + yTran;
% 	haloPos(:,3) = haloPos(:,3) + zTran;
% 
% 	% Scale
% 	scaleX = transformParams(7);
% 	scaleY = transformParams(8);
% 	scaleZ = transformParams(9);
% 
% 	% Rotate
%     thetaX = transformParams(4);
%     thetaY = transformParams(5);
%     thetaZ = transformParams(6);
%     Rx = [1, 0, 0; 0, cos(thetaX), -sin(thetaX); 0, sin(thetaX), cos(thetaX)];
%     Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
%     Rz = [cos(thetaZ), -sin(thetaZ), 0; sin(thetaZ), cos(thetaZ), 0; 0, 0, 1];
% 
% 	
% 
% 	% Second alignment
% 	C = cadPos' * haloPos;
% 	[U, ~, V] = svd(C);
% 	R = V * U';
% 	if det(R) < 0
% 		V(:,3) = -V(:,3);
% 		R = V * U';
% 	end
%     Ropt = Rx * Ry * Rz * R;
%     transformedHaloPos = haloPos * Ropt;
% 
% 	% Compute error
%     distances = sqrt(sum((transformedHaloPos - cadPos).^2, 2));
%     error = mean(distances);
% end

%==========================================================================
% - G E T   C H A N N E L   N U M B E R   F R O M   L V M
%==========================================================================
function positions = getchannelNumberFromLvm(lvm,positions)
%-READ channel information from lvm file and rename positions
%--------------------------------------------------------------------------
% Code is adapted from info2pos_neuro1 by Stephanie Mellor

% Open lvm
[~, lvmFile] = fileparts(lvm);
disp(['Opening ', lvmFile])
fid = fopen(lvm);
for lin = 1:18
	line = fgetl(fid);
end
units = textscan(line, '%s', 'Delimiter', '\t');
units = units{1}(2:end);
if any(cellfun(@isempty, units))
	idx = cellfun(@isempty, units);
	units(idx) = {'other'};
end
for lin = 19:23
	line = fgetl(fid);
end
channels = textscan(line, '%s', 'Delimiter', '\t');
channels = channels{1}(2:end);                  % Cut out time channel
channels(startsWith(channels, 'Comment')) = []; % Remove empty comment channel
channels(startsWith(channels, 'Untitled')) = {'DataLogger'};
fclose(fid);

% Find channels overlapping with slot2sens.csv
meg_chans = contains(units, 'pT');
channels = channels(meg_chans);

% Get rid of duplicates (X,Y,Z)
if ~any(contains(channels, '_'))
	% This assumes the .lvm channels are labeled 'X9' etc. and are all G3
	numbered_naming = true;
elseif all(cellfun(@length, channels) == 4)
	% This assumes the .lvm channels are labeled 'B1_X' etc. and are all G3
	numbered_naming = false;
else
	% This assumes the .lvm channels are labeled 'B9_X' etc. and are all G3
	% Relabel to have format X9 etc.
	numbered_naming = true;
end

% Relabel channel names in positions from B1 to 9, if names in lvm file are labelled as such
if isa(positions.name, "cell") && numbered_naming
	board_slots = cell(size(positions, 1), 1);
	for slot = 1:size(positions,1)
    	if isa(positions{slot, 'name'}{1}, "char")
        	board_slots{slot} = regexp(positions{slot, 'name'}{1}, '[A-H]', 'match');
        	if isa(board_slots{slot}, "cell")
            	board_slots{slot} = board_slots{slot}{1};
        	end
    	end
	end
	rename_inds = cellfun(@isempty, board_slots);
	rename_inds = ~rename_inds;
	board_numbers = zeros(size(board_slots));
	board_numbers(rename_inds) = cellfun(@(x)double(upper(x))-double('A')+1, board_slots(rename_inds));
	chans = zeros(size(board_numbers));
	chans(rename_inds) = (board_numbers(rename_inds)-1)*8 + cellfun(@(x)str2double(x(2)), positions.name(rename_inds));
	chans(~rename_inds) = cellfun(@(x)str2double(x), positions.name(~rename_inds));
end

% Rename the positions table
for posIdx = 1:length(positions.name)
	curChan = num2str(chans(posIdx));
	positions(posIdx,:).name{1} = [curChan,'-',positions(posIdx,:).name{1}];
end
end

%==========================================================================
% - G E T   C H A N N E L   N U M B E R   F R O M   L V M
%==========================================================================
function positionNumber = getChannelNumberFromLvm(boardNameNumber)
letter = boardNameNumber(1);
number = str2double(boardNameNumber(2));
letterIndex = double(upper(letter)) - double('A');
numberIndex = number - 1;
positionNumber = bitshift(letterIndex, 3) + numberIndex + 1;
end