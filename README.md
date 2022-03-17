## Repository for the Tomogram simulation based on PyTom

### Done.
1. volume object(and .em file) generator from given PDB ID(.pdb, .cif both supported).   
2. compact volume generator(Cuboidal form).   
3. mrc file formatter
4. naive multiparticle generator.   
   1. Crowd scenario generator
   2. Particle incorporator
5. test rotation scheme   
   1. Compact Cuboid
   2. Relaxation to Minimum Cube
   3. Rotate
   4. Re-compact Cuboid
6. random rotation with middle efficiency.
7. PYTOM's em2mrc generates mrc file much faster. [ for certain file(1_withoutrotation.em) : 99s vs 6.6s ]
   1. But, mrcfile.validate -> false.
   2. int mode mrc generated. -> File size is efficient : 134.2M vs 536.9M
8. Dimension adjustment
   1. maybe just 1 A -> 10 A(1 nm)

### Working.   
1. Rotation : with stride parameter [ saved ] / with stride parameter [ unsaved ]
2. Image formation options...
3. multi particle generator?
   1. crowd scenario generator   
   2. scenario to volume file converter
4. Other processing options?

### Timeline
20220227 00:50 - compact volume generator completed.      
20220306 23:47 - Update Readme for rotation TODO.   
20220313 24:00 - (Until today) Rotation, Resolution adjustment done.   
20220314 
