with import <nixpkgs> {};

stdenv.mkDerivation rec {
  name = "env";
  env = buildEnv { name = name;
                  paths = buildInputs; };
  buildInputs = [
    gcc6
    fftw
    hdf5
    valgrind
    fftwFloat
    python35
  ];
}


