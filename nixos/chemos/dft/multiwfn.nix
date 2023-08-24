{ stdenv
, lib
, autoPatchelfHook
, unzip
, fetchurl
, lapack
}:

stdenv.mkDerivation rec {
  pname = "multiwfn";
  version = "3.8";

  src = fetchurl {
    url = "http://sobereva.com/multiwfn/misc/Multiwfn_3.8_dev_bin_Linux_noGUI.zip";
    sha256 = "sha256-Z+JnpUreVX0imKmRoUsLRBQDtog8xE5n/0AmJVXpf+c=";
  };

  nativeBuildInputs = [
    autoPatchelfHook
    unzip
  ];

  buildInputs = [
    lapack
  ];

  sourceRoot = ".";

  installPhase = ''
    install -m755 -D Multiwfn_${version}_dev_bin_Linux_noGUI/Multiwfn_noGUI $out/bin/multiwfn
  '';

  meta = with lib; {
    homepage = "http://sobereva.com/multiwfn/";
    description = ''
      Multiwfn is an extremely powerful program for realizing electronic
      wavefunction analysis, which is a key ingredient of quantum chemistry.
    '';
    platforms = platforms.linux;
    mainProgram = "multiwfn";
  };
}
