with import <nixpkgs> {};

stdenv.mkDerivation rec {
  name = "env";
  env = buildEnv { name = name;
                  paths = buildInputs; };
  buildInputs = [
    gcc-snapshot
    fftw
    hdf5
    astropy
    graphviz
    kdeApplications.kcachegrind
    python35
    python35Packages.pylru
    python35Packages.h5py
    python35Packages.scipy
    python35Packages.numpy
    valgrind
    fftwFloat
  ];
}


