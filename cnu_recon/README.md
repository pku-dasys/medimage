# 3D Industrial CT Reconstruction

# test datasets

We now only support test datasets using parallel beam.

1. Run `./bin/generate` to generate `test.dr`.
2. Run `./bin/ct3d test.json` to write image results in `test` folder.
3. You can use Octave or Matlab to show the images.
4. Or run `./viewer/viewer test_dir NX NY NZ` to show the images.

# issue
1. Only works well in model 1 & 2, but not in model 3 (4-layer ellipsoid)

