import numpy as np

def overlay_qmri_over_anatomical(img,mask,map_vals,map_min,map_max,colorinputs,bightness_val=1.0):
	''' Creates a RGB image where qMRI is overlaid onto anatomical MRI
	qMRI_over_anat = ...
           ... overlay_qmri_over_anatomical(img,mask,map_vals,map_min,map_max,colorinputs,brightness_val=1.0)
	
	INPUTS
	 1) img: anatomical image stored as 2D numpy floating point array 
	          (it can be in any range and is converted to
	          numpy unit8, s.t. the minimum vale of img correponds to 0 and the
	          maximum value to 255)
	 2) mask: a binary mask having the same size of img stored as a 2D numpy array.
	          It should be set to 0 for pixels 
	          that in the output image will show information from 
	          img, 1 for those that will show information for parametric map
	 3) map_vals: a 2D numpy array storing the qMRI map to be overlaid onto the 
	              anatomical information, having the same size as img.
	          error
	 4) map_min: the value of map to be mapped to the first colour of
	              colour list stored in colorinputs. NaNs are mapped to map_min
	 5) map_min: the value of map to be mapped to the last colour of
	              colour list stored in colorinputs
	 6) colorinputs: a 2D numpy matrix of Ncolours x 3 elements storing the colour list
	                 to which values in map will be mapped linearly (first
	                 coloumn: red; second column: blue; third column: green.
	                 The RGB coordinate can be expressed in [0;1]x[0;1]x[0;1]
	                 or [0;255]x[0;255]x[0;255] as integers)
	 7) bightness_val: numerical value to which the maximum value in img will
	                   be mapped (bightness_val > 0); default 1.0, set it to a
	                   number < 1.0 to get a darker anatomical image or to a 
	                   number > 1.0 to get a brighter anatomical image
	 
	 OUTPUTS
	 1) qMRI_over_anat: 3channel RGB image (3D uint8 numpy array)
	
	 ALGORITHM: values of the qMRI metric are mapped linearly to the cololur
	 list stored in colourinputs.
	
	Author: Francesco Grussu, University College London
		   <f.grussu@ucl.ac.uk> <francegrussu@gmail.com>'''

	# Convert anatomical image to uint8
	anat = np.uint8( 255*(bightness_val*img - np.min(img)) / (np.max(img) - np.min(img)) )
	
	# Convert 2D floating point matrix of colours to a 2D uint8 matrix	
	colorset = np.uint8(255*colorinputs)
	Ncolours = colorset.shape[0]
	
	# Pointer mapping linearly the value of the qMRI metric to the set of colours available
	metric_index = np.linspace(map_min,map_max,Ncolours)  # Value of the qMRI metric corresponding to each colour
	
	# Create the output image: initialise it as anatomical
	qMRI_over_anat = np.dstack((anat,anat,anat));
	rows = qMRI_over_anat.shape[0]
	cols = qMRI_over_anat.shape[1]
	
	# Loop over pixels to put colours for the qMRI metric
	for rr in range(0,rows):
		for cc in range(0,cols):
		
			# Add parametric map with colours if within the mask
			if(mask[rr,cc]):
				
				# Get value of qMRI metric at pixel (rr,cc)
				qMRIval = map_vals[rr,cc]
				if(np.isnan(qMRIval)):    # if metric is NaN, replace it with the minimum value of the colorbar
					qMRIval = map_min

				# Find the metric value most similar to it
				colourDiff = np.abs(metric_index - qMRIval)
				colourPointer = np.argmin(colourDiff)
				       
				# Get the corresponding colour
				qMRIcolour = colorset[colourPointer,:]
				
				# Save the colour in the image to be plotted
				qMRI_over_anat[rr,cc,0] = qMRIcolour[0]    # Red channel
				qMRI_over_anat[rr,cc,1] = qMRIcolour[1]    # Green channel
				qMRI_over_anat[rr,cc,2] = qMRIcolour[2]    # Blue channel
	
	# Return output uint8 RGB image with parametric map in colour overlaid onto grey-scale anatomical scan
	return qMRI_over_anat
	

