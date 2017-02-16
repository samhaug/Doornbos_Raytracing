# Doornbos_Raytracing

The command is ./ray1d

You will have to change line 31 of modps.f. Set this line to the directory in which rmod.dat lives. 
The code needsto find this file. 

rm *.o and make

The input file is ray1d.in. Run the executable like this:

./ray1d < ray1d.in

Sample entry explanation for SAMPLE.in

SAMPLE.in has many ‘blocks’, but the code will only run on the first block in the input file. You may specify the number of phases by the number in the first line.

1 is number of phases computed, 0.5 is a stabilization number that is less than 1. Satoshi has never changed it.
1   0.5 
5.2433 is ballpark guess for slowness, 2 is for lines that follow                                                                  
5.2433          2
A simple phase like P has two branches, the downgoing part and the upgoing part. Hence two lines:
5871.0 is radius for source, 3479.5 is radius for turning point. 
For P wave, make this the CMB. The middle 1 means P.                                                             
5871.0  3479.500000  1 1  1
This is the upgoing part, modeled as a source at the earth’s surface and a receiver at the CMB, so it goes backwards, hence the -1                                                    
6371.0  3479.500000  1 1 -1                                                    
 -- P for h=500 D=80  AK135

This is a similar block of text for PcP. Everything is the same, except the slowness guess on the second line. This dictates what raypath the code searches for.
     1   0.5                                                                        
     4.37606         2                                                              
     5871.0  3479.500000  1 1  1                                                    
     6371.0  3479.500000  1 1 -1                                                    
      -- PcP for h=500 D=80  AK135





This block is for P4KP. There are three branches to specify. The in-going P, the 4K, and the outgoing P.
      1   0.5                                                                        
      4.42815  3                 
The in-going P is the same as usual, it has one branch                                                 
      6371.0  3479.500000  1 1 1
The core bouncing branch has 4 legs, each leg has a downward and upward going section, making 8 total. Put the receiver at the center of the earth.                                                                                         
      3479.5  0.000000  8 1 1     
      The out-going P is the same as usual, it has one branch going backwards                                                
      6371.0  3479.500000  1 1 -1                                                    
      -------- P4KP ----------  

ScS PREM block. Top two lines are the same as previous parts.
      1   0.5                                                                        
      6.93724966      2
First branch, going backwards. Starting at surface and going to CMB. The middle number is 2 to indicate an S phase.                                                                     
      6371.0  3480.000000  1 2 -1      
Second branch going forward.                                
      6171.0  3480.000000  1 2  1                                                          
      ------ ScS for PREM   ----

P to S conversion. The P branch only goes downward, where it converts to a S wave. The S wave has a downgoing and upgoing part. So three lines to describe the entire phase.
      1   0.5                                                                        
      9.80711124       3                                                             
      5731.0 5711.0  1  1  1                                                        
      5711.0 3482.0  1  2  1                                                        
      6371.0 3482.0  1  2 -1                                                        
      ----- P660S-h640km 84 deg-----

      Doing three
      3   0.5                                                                                                            
      3.5                2                                                             
      5741.0  3482.000000  1 1  1                                                    
      6371.0  3482.000000  1 1 -1                                                    
      3.331            3                                                             
      6371.0  5741.000000  1 1 -1                                                    
      6371.0  3482.000000  1 1  1                                                    
      6371.0  3482.000000  1 1 -1                                                    
      3.372            3                                                             
      6371.0  5741.000000  1 2 -1                                                    
      6371.0  3482.000000  1 1  1                                                    
      6371.0  3482.000000  1 1 -1                                                    
      ---PcP, pPcP, sPcP for IASP91 ----
