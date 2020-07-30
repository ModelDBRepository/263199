This is the readme for the models associated with the paper:

Smith P, Buhl E, Tsaneva-Atanasova K, Hodge JJL (2019)
Shaw and Shal voltage-gated potassium channels mediate circadian changes in 
Drosophila clock neuron excitability. J Physiol



LNVmodel models the activity of a Drosophila LNV neuron
Necessary files: LNVmodel	(main code)
		 MODEL.mat	(model parameters)

		 



To run the code on Matlab:

Call the function:
LNVmodel(ZT, kv2choice);

where ZT is the time of day  i.e. ZT0 is lights-on, ZT12 is lights-off 
(in a 12hr light/12hr dark cycle), ZT6 is mid-day etc.

where kv2choice is the choice of kv2 models: 
	1 is native kv2
	2 is native kv2 with human wild-type kv9 also expressed
	3 is native kv2 with human mutant kv9 (c379E) also expressed





To return outputs use:
[vrec, I1, I2, I3, I4] = LNVmodel(ZT, kv2choice);

Where vrec is the membrane voltage of the modelled neuron
      I1 is the modelled current of the Shaker (Kv1) channel
      I2 is the modelled current of the Shab (Kv2) channel
      I1 is the modelled current of the Shaw (Kv3) channel
      I1 is the modelled current of the Shal (Kv4) channel





To create different figures:

The default output is a graph (figure 5) of the membrane voltage of the 
modelled neuron (vrec) along with a graph of the currents of each of the 
four channels Shaker (I1), Shab (I2), Shaw (I3), and Shal (I4). To reproduce
figure 5 of the paper run the model with kv2choice=1 and ZT=0 (morning) or 
ZT=12 (evening), changing subplots and x-axis as required. 

To change figures go to the 'Creating figures' subsection towards the end 
of the script (Line 223).

