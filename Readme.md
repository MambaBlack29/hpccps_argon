
**Compile**

type  'make ' in this folder to compile the codes. On successful installation, it creates the executable 'mdrun' 

**clean the folder for re-installation**
type  'make clean' 

**How to run**
1. MD input file  'md.input'
2. Modify the input file 'md.input' to change  a) no of molecules b) input filename c) box length 
	For liquid argon, use the following:
		Atoms		Box length(A)
		1600		42.3212
		5400		63.4818
		12800		84.6420

3. To run the simulation, type the following: './mdrun md.input' 

4. To know the runtime, check the last line of md.out file 

