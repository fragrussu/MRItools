function Temp = WaterADC2Temperature(ADC)
% Converts and ADC map of water into a temperature map 
%
% Temp = WaterADC2Temperature(ADC)
%
% INPUT: ADC (map of ADC of water, expressed in m^2 / sec).
% OUTPUT: Temp (temperature, expressed in Celsius).
%
% The function calculates the temperature T from self-diffusivity D as
%
%        T =  215.05*(1 + (D/1.635e-8)^(1/2.063)) - 273.15;
%
% as shown in Holz et al, "Temperature-dependent self-diffusion coefficients 
% of water and six selected molecular liquids for calibration in accurate 
% 1H NMR PFG measurements", Phys Chem Chem Phys (2000), vol 2, pag 4740-4742.
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

Temp = 215.05*(1 + (ADC/1.635e-8).^(1/2.063)) - 273.15;

end
