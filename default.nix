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
    python35Packages.h5py
    python35Packages.scipy
    python35Packages.numpy
    python35Packages.pylru
    valgrind
    fftwFloat
    python35
  ];
}


