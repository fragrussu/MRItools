function ReadTrajFilesCamino(trajfile,datafile)
% Reads a trajectory file from the datasynth command of Camino for PGSE
% 
% ReadTrajFilesCamino(trajfile,datafile)
%
% INPUT PARAMETERS 
% trajectory  -->   trajectory file from Camino
% datafile    -->   output data file
%
% datafile is a .mat file containing the variable data, which is a 
% structure with following fields:
% 
%         data.TE --> echo time
%         data.dt --> time step
%         data.Nspins --> number of spins
%         data.Nsteps --> number of time steps, equal to TE/dt
%         data.X  --> x position of spins (matrix of size NstepsxNspins)
%         data.Y  --> y position of spins (matrix of size NstepsxNspins)
%         data.Z  --> z position of spins (matrix of size NstepsxNspins)
%         data.T --> time (array of size Nsteps x 1)
%
% If saving a structure is not allowed in the Matlab version employed to
% run the function, then the different fields are saved as separate
% variables.
%
% Author: Francesco Grussu, University College London <f.grussu@ucl.ac.uk>
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

% Open file for reading
fid = fopen(trajfile,'r','b');   % Open file for reading as big-endian

fprintf('Reading file: %s\n',trajfile)

% Read simulation duration
TE = fread(fid,1,'double');

% Read number of spins
Nspins = fread(fid,1,'double');

% Read number of time steps
Nsteps = fread(fid,1,'double');


% Get duration of a single step
dt = TE/Nsteps;

% Allocate trajectories
X = zeros(Nsteps,Nspins);
Y = zeros(Nsteps,Nspins);
Z = zeros(Nsteps,Nspins);
T = zeros(Nsteps,Nspins); 
I = zeros(Nsteps,Nspins);

% Read data
for tt=1:Nsteps
    
    
    
    fprintf('          time step %d out of %d. Please wait, reading %d spins...\n',tt,Nsteps,Nspins);
    for ss=1:Nspins
        
        T(tt,ss) = fread(fid,1,'double');
        I(tt,ss) = double(fread(fid,1,'int'));   % Spin index for current datum
        X(tt,ss) = fread(fid,1,'double');
        Y(tt,ss) = fread(fid,1,'double');
        Z(tt,ss) = fread(fid,1,'double');
        
    end
    fprintf('\n');
end
clear I
fclose(fid);

fprintf('\n\nData read. Now saving to disk...\n\n')

% Save data to disk
data.TE = TE;
data.dt = dt;
data.Nspins = Nspins;
data.Nsteps = Nsteps;
data.X = single(X); clear X
data.Y = single(Y); clear Y
data.Z = single(Z); clear Z
data.T = single(T(:,1)); clear T

try 
    save(datafile,'data','-v7.3');
catch
    X = data.X;
    Y = data.Y;
    Z = data.Z;
    T = data.T;
    TE = data.TE;
    dt = data.dt;
    Nspins = data.Nspins;
    Nsteps = data.Nsteps;
    clear data    
    save(datafile,'X','Y','Z','T','TE','dt','Nspins','Nsteps','-v7.3');     
end

fprintf('\n\nDone...\n\n')

end



