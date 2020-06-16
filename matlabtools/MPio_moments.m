function MPio_moments(input,output,noise,kernel,varargin)
% Run MPdenoising() on an NIFTI file and mitigate noise floor
% 
% MPio_moments(input,output,noise,kernel)
% MPio_moments(input,output,noise,kernel,mask)
%
% input:   path to a NIFTI file storing the scan to denoise
% output:  path of output NIFTI file storing the denoised scan
% noise:   path of output NIFTI file storing the noise map 
% kernel:  3x1 (or 1x3) array storing the kernel size for 
%          the denoising sliding window
% mask:    path of an input NIFTI file with a mask excluding
%          areas where denoising is not required
%
%
% Dependencies:
% MPdenoising() Matlab function by Veraart J et al
% https://github.com/NYU-DiffusionMRI/mppca_denoise/blob/master/MPdenoising.m
%
%
% References
%
% - MP-PCA ALGORITHM:
% Veraart J et al, Neuroimage 2016, 142:394-406 
% DOI: 10.1016/j.neuroimage.2016.08.016
%
%
% - NOISE FLOOR MITIGATION:
% Koay CG and Basser PJ, J Magn Reson 2006, 179(2):317-22
% DOI: 10.1016/j.jmr.2006.01.016 
%
% Author: Francesco Grussu <f.grussu@ucl.ac.uk>
%
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

    % Load input NIFTI
    fprintf('    ... loading input data\n');
    info = niftiinfo(input);
    hin = niftiread(info);
    hin = double(hin);
    
    % Denoise with MP-PCA
    fprintf('    ... denoising with MP-PCA\n');
    if(isempty(varargin))
        [hden,sgmmap] = MPdenoising(hin, [], kernel);
    else
        tissuemask = niftiread(varargin{1});
        tissuemask = logical(tissuemask);
        [hden,sgmmap] = MPdenoising(hin, tissuemask, kernel);        
    end
    

    % Mitigate noise floor with method of moments
    fprintf('    ... mitigating noise floor with the method of moments\n');
    hden_corr = zeros(size(hden));
    sgmmap_corr = zeros(size(hden));

    for mm=1:size(hden,1)
        for nn=1:size(hden,2)
            for kk=1:size(hden,3)
                
                if(~isempty(varargin))
                    maskval = tissuemask(mm,nn,kk);                   
                else
                    maskval = true;
                end
                
                if(maskval)
                    for qq=1:size(hden,4)
                        % Mitigate noise floor in each measurement
                        data_den_meas = hden(mm,nn,kk,qq);
                        s_r = sgmmap(mm,nn,kk);
                        snr_opt = iteratesnrmap(data_den_meas,s_r,1e-08);
                        corr_opt = snrfunc(snr_opt);
                        sgmmap_corr(mm,nn,kk,qq) = sqrt(s_r*s_r/corr_opt);
                        data_den_meas_new = sqrt( data_den_meas*data_den_meas + s_r*s_r*(corr_opt - 2)/corr_opt );   % Mitigate floor
                        if(isreal(data_den_meas_new))   % Check that signal did actually show presence of noise floor
                           hden_corr(mm,nn,kk,qq) = data_den_meas_new;
                        else
                            hden_corr(mm,nn,kk,qq) = data_den_meas; % Measure lower than noise floor (compatible with Guass. noise) 
                        end
                    end
                end
            end            
        end
    end

    % Get one summary value for sigma of noise in each voxel
    sgmmap_corr_mean = mean(sgmmap_corr,4);

    % Save denoised NIFTI
    fprintf('    ... saving output data\n');
    infoout = info;
    infoout.Datatype = 'double';
    niftiwrite(hden_corr, output, infoout);
    
    % Save noise map
    infonoise = info;
    infonoise.Datatype = 'double';
    infonoise.ImageSize = info.ImageSize(1:3);
    infonoise.PixelDimensions = info.PixelDimensions(1:3);
    niftiwrite(sgmmap_corr_mean, noise, infonoise);

end







function [theta_opt,niter] = iteratesnrmap(mag_expect,mag_std,tol)
%
% Discrete map to be used for noise floor correction
%
% mag_expect: expected value of magnitude measurements
% mag_std: standard deviation of magnitude measurements
%
% Please check DOI: 10.1016/j.jmr.2006.01.016
%

	% Maximum number of iterations
	NITERMAX = 80;
	% Iteration counter
	niter = 1;
	% Lower bound for the ratio of moments
	minval = sqrt(pi/(4.0-pi));
	% Get ratio of moments
	momratio = mag_expect/mag_std;

	% Check whether the ratio is plausible or not
	if(momratio<=minval)
	    theta_opt = 0;
	else

		% Initialise the iterative map
		theta = momratio - minval;	
		g_of_theta = snrmap(theta,mag_expect,mag_std);

	       % Iterate map
		while niter<NITERMAX
			tolval = abs(theta - g_of_theta);
			if tolval<tol
			    theta_opt = theta;
			    break
			else
				theta = g_of_theta;
				g_of_theta = snrmap(theta,mag_expect,mag_std);
				theta_opt = g_of_theta;
				niter = niter + 1;
			end
	    	end
	    
	end


end






function eps_of_theta = snrfunc(theta)
% Function evaluating a correction factor as a function of SNR
% Please check DOI: 10.1016/j.jmr.2006.01.016
%

	    if(theta<35)
		eps_of_theta = 2 + theta.*theta - (pi/8)*exp(-0.5*theta.*theta).*( (2 + theta.*theta).*besseli(0,0.25*theta.*theta) + theta.*theta.*besseli(1,0.25*theta.*theta)).*( (2 + theta.*theta).*besseli(0,0.25*theta.*theta) + theta.*theta.*besseli(1,0.25*theta.*theta));
	    else
		eps_of_theta = 1;
	    end
		

end




function g = snrmap(theta,mag_expect,mag_std)
% Function evaluating an iterative map used for SNR calculation
%
% theta: SNR (true signal / gaussian noise std)
% m_mag: expected value of rician-biased measurements
% m_std: std of rician-biased measurements
%
% Please check DOI: 10.1016/j.jmr.2006.01.016
%

	g = sqrt(snrfunc(theta)*(1 + mag_expect*mag_expect/(mag_std*mag_std)) - 2);


end

