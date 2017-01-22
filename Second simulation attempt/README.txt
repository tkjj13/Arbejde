Main files
	Simulation_v1_0
		Possible to change number of carriers but with same BW much like LTE still very crude and with a few critical flaws
	Simulation_v1_1
		Corrected most of the critical flaws but still same BW
	Simulation_v1_2
		First attempt at getting different BW for the carriers
		
		
		
		
		
Function files
	ser2par
		Takes vector input and split to a matrix
	symbolGen
		Generate symbols from bit sequence
	symbolDegen
		Generate bit from symbol sequence using hard ML detect
		
		
Unused files
	myBin2dec
		Takes binary vector convert to decimal number
	myDec2bin
		Takes decimal number convert to binary vector
	RFchannel
		Takes RF input and simulate channel (NOT COMPLETE)
	rfDegen
		Down convert RF to CBB (NOT COMPLETE)
	rfgen
		Up convert CBB to RF (NOT COMPLETE)
	coding
		Encode bit stream for crypto, source and channel (NOT COMPLETE)
	Bit_fil_generator
		Generate random test bits and place in Data.MAT 
	SymbolChannel
		Takes symbol input and simulate channel (NOT COMPLETE)


.MAT files
	Data
		Contains random binary vectors