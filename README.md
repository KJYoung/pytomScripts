## Repository for the Tomogram simulation based on PyTom
This scripts are under development. Not optimized code.

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
6. Random rotation with rotationStep.
7. PYTOM's em2mrc generates mrc file much faster. [ for certain file(1_withoutrotation.em) : 99s vs 6.6s ]
   1. But, mrcfile.validate -> false.
   2. int mode mrc generated. -> File size is efficient : 134.2M vs 536.9M
8. Dimension adjustment
   1. maybe just 1 A -> 10 A(1 nm)?
   2. TODO : re-verify.
9. Subtomogram generator


### TODO.   
1. Image formation options...
2. Minus Noise
3. Subtomo include?
4. Wrong dimension adjustment?
5. Other processing options?

### Timeline
20220227 00:50 - compact volume generator completed.      
20220306 23:47 - Update Readme for rotation TODO.   
20220313 24:00 - (Until today) Rotation, Dimension adjustment done.   
20220321 16:00 - (Until today) Subtomogram generator, Metadata writer done.
