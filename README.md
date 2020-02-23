# US2 to SE File Converter
Quick Description: Converts .simulation files from Universe Sandbox 2 to .sc file for use with Space Engine 0.99



Universe Sandbox 2 to Space Engine 0.99
			File Converter 1.0
				Created by Seth Fusion
				
				
Quick Start Instructions:

	1. For a video walkthrough, go here: (Watch the whole video)
		https://youtu.be/0Z4PdGNjLMI
	
	2. Find the simulation file you want to convert from your US2 documents folder.
		This is probably found in your documents folder, but the easiest way to get
		to it is to open US2, click on "Open Existing Simulation", click on "My Sims", 
		and click "Open Simulation Folder". This folder contains .ubox files, which 
		can be opened with 7zip or anything like it. Open the .ubox folder for the
		simulation you want to convert.
		
	3. Inside the .ubox folder, the only file you will need is the .simulation file. 
		Copy/paste that file into the "input" folder that this file came with. 
		
	4. Come back out of the "input" folder and double click on the "US2 to SE Converter.exe"
		A black box will appear, but you don't have to do anything. Depending on the size of 
		your simulation, this may take a few seconds, but not too long. 
		
	5. When the black box goes away, open the "output" folder and find the two files that
		were just created. One file will end with the name "planet.sc", and the other
		will end with "star.sc". The star file goes under the "stars" folder of your SE addons.
		The planet file goes under your planet folder in you SE addons. 
		
	6. For the general case, that should be it. You can open Space Engine and search for your
		system!
		
	If there are any issues, other than the disclaimer below, please message me on the forums
		and I will look into fixing it. 
		
IMPORTANT DISCLAIMER please read this

	Binary orbits DO NOT work. At least, not completely. The program will generally calculate
	the semi-major axis correctly (mostly), but all other values are borderline useless. 
	I wish there was something I could do about it, but I'm not sure if the problem is with 
	my code, or if binary orbits are too complex to calculate with a single instance of information. 
	
	Even though they don't work, SE will usually still render the objects, and it is still a good
	staring point. But, you will need to open the .sc file, find the objects, and fix their orbits 
	manually if you want them to look normal. You can do this by fixing the values for one object 
	and copying/pasting the values over to the second object, but remember to subtract 180 from ONE 
	object's ArgOfPericenter. 
