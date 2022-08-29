# Intrinsic-Cyclizability-Prediction-Codes
Codes to predict intrinsic cyclizability of DNA sequences
These codes are associated with the publication "Deciphering the mechanical code of the genome and epigenome", and are referred to in the supplementary notes.

Predict_linear.m: MATLAB code to predict intrinsic cyclizability of sequences on the basis of the model where intrinsic cyclizability is a linear combination of the 16 dinucleotide contents.

Predict_periodic.m: MATLAB code to predict intrinsic cyclizability of sequences on the basis of the model where intrinsic cyclizability is a linear combination of the 136 Helical Repetition Indices (HRIs).

Predict_combined.m: MATLAB code to predict intrinsic cyclizability of sequences on the basis of the model where intrinsic cyclizability is a linear combination of the 16 dinucleotide contents and the 136 HRIs.

NN_sequence.m: MATLAB code that contains the function "NN_sequence", which takes as imput a 200 element vector, which represents a 50 bp sequences via one hot coding, and produces as output the predicted intrinsic cyclizability, based on a neural nets model.

NN_152param.m: MATLAB code that contains the function "NN_151param", which takes as input a vector of 151 elements (15 independent dinucleotide contents and 136 HRIs), which represents a 50 bp DNA sequences, and produces as output the predicted intrinsic cyclizability, based on a neural nets model.
