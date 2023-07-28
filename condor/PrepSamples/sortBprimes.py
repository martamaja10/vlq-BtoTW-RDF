import os,sys,shutil,datetime,time

sample_string = ""

with open('BPrimeNanos.txt', 'r') as file:
    # Method 1: Using readline()
    line = file.readline()
    while line:
        # Process the line
        mass = line[15:19]
        
        year = ""
        if "20UL16" in line:
            year = "20UL16"
        elif "20UL17" in line:
            year = "20UL17"
        elif "20UL18" in line:
            year = "20UL18"
        
        inputFile = "Bprime_M" + mass + "_" + year
        if "20UL18" in line:
            sample_string = sample_string + "    \"" + inputFile + "\":" + inputFile + ",\n"
                
        #line = line.strip()
        #sample_string = sample_string + inputFile + " = sample(\"" + inputFile + "\", \"" + inputFile + "NanoList.txt\", \"" + line + "\")\n"

        # Read the next line
        line = file.readline()
        
with open('BPrimeSamples.txt', 'w') as file:
    file.write(sample_string)