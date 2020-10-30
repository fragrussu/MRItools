function qMRI_over_anat = overlay_qMRI_over_anatomical(img,mask,map,map_min,map_max,colorinputs,bightness_val)
% Creates a RGB image where qMRI is overlaid onto anatomical MRI
%
% qMRI_over_anat = ...
% ... overlay_qMRI_over_anatomical(img,mask,map,map_min,map_max,colorinputs,brightness_val)
%
% INPUTS
% 1) img: anatomical image (it can be in any range and is converted to
%          unit8(), s.t. the minimum vale of img correponds to 0 and the
%          maximum value to bightness_val; it should NOT be an integer image)
% 2) mask: 0 for pixels that in the output image will show information from 
%          img, 1 for those that will show information for map
%          (a binary mask having the same size of img)
% 3) map: the qMRI map to be overlaid onto the anatomical information
%          (a 2D matrix of the same size of anat) -- NaN values will give
%          error
% 4) map_min: the value of map to be mapped to the first colour of
%              colour list stored in colorinputs
% 5) map_min: the value of map to be mapped to the last colour of
%              colour list stored in colorinputs
% 6) colorinputs: a matrix of Ncolours x 3 elements storing the colour list
%                 to which values in map will be mapped linearly (first
%                 coloumn: red; second column: blue; third column: green.
%                 The RGB coordinate can be expressed in [0;1]x[0;1]x[0;1]
%                 or [0;255]x[0;255]x[0;255] as integers.
%                 For instance, one could use 
%                           colorinputs = colormap('jet');      )
% 7) bightness_val: numerical value to which the maximum value in img will
%                   be mapped (bightness_val > 0); suggested option 255
% 
% OUTPUTS
% 1) qMRI_over_anat: 3channel RGB image to be shown with imshow or saved
%                    with imwrite
%
% ALGORITHM: values of the qMRI metric are mapped linearly to the cololur
% list stored in colourinputs.
%
% Author: Francesco Grussu, University College London 
%         <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>
%
% BSD 2-Clause License
% 
% Copyright (c) 2020 University College London
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 

% Check that the brightness value makes sense
if((bightness_val<0)||(isscalar(bightness_val)==0)||(isnumeric(bightness_val)==0))
    error('overlay_qMRI_over_anatomical(): the brighness input bightness_val should be a scalar, numerical value greater than 0.') 
end

% Convert anatomical image to single-channel gray scale (8-bit integers)
if(isinteger(img))
   error('overlay_qMRI_over_anatomical(): input image img must be a single/double-precision image and not an integer image.') 
end
anat = uint8( bightness_val*(img - min(min(img)) ) / ( max(max(img)) - min(min(img)) ) );

if(isinteger(colorinputs)==0)
    colorset = uint8(255*colorinputs);
else colorset = colorinputs;
end

% Number of colours available to plot the qMRI map
Ncolours = size(colorset,1);

% Pointer mapping linearly the value of the qMRI metric to the set of colours available
metric_index = linspace(map_min,map_max,Ncolours);  % Value of the qMRI metric corresponding to each colour

% Create the output image: initialise it as anatomical
qMRI_over_anat = cat(3,anat,anat,anat);
[rows,cols,~] = size(qMRI_over_anat);

% Loop over pixels to put colours for the qMRI metric
for rr=1:rows
    for cc=1:cols
       
        % check whether the pixel is within the mask...
        if(mask(rr,cc)==1)   % ... yes, it is. We need to find the best colour representing the metric at this pixel
           
            % Get value of qMRI metric at pixel (rr,cc)
            qMRIval = map(rr,cc);
            if(isnan(qMRIval))    % if metric is NaN, replace it with the minimum value of the colorbar
                qMRIval = map_min;
            end
            % Find the metric value most similar to it
            colourDiff = abs(metric_index - qMRIval);
            colourPointer = find(colourDiff==min(colourDiff));
            colourPointer = colourPointer(1);
       
            % Get the corresponding colour
            qMRIcolour = colorset(colourPointer,:);
            % Save the colour in the image to be plotted
            qMRI_over_anat(rr,cc,1) = qMRIcolour(1); % Red channel
            qMRI_over_anat(rr,cc,2) = qMRIcolour(2); % Green channel
            qMRI_over_anat(rr,cc,3) = qMRIcolour(3); % Blue channel
        end
        
    end 
end


end
