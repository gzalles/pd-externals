pool_vox~
======
this is my external for the pool instrument 

it takes in 3 floats (HSV) values and returns 16 ambisonic channels of encoded audio (3OA) based on creation args azi and elev. 

the hue is used to determine the frequency of each voice, saturation determines the cutoff freq of a low-pass filter, and value is used to change the overall volume. specifically (1-value) controls the amplitude of each voice.  

