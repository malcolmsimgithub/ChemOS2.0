{ buildPythonPackage
, fetchPypi
, lib
, numpy
, scipy
}:

buildPythonPackage rec {
  pname = "pyDOE";
  version = "0.3.7";
  format = "setuptools";
  doCheck=false;

  src = fetchPypi {
    inherit pname version format;
    sha256 = "fac921ef31e4642f86c2e5ea2b1fe8c3d267b44259171eb4be9c704be6a7f489";
  };

  propagatedBuildInputs = [
    numpy
    scipy
  ];

  meta = with lib; {
      description = "Design of experiments for python";
      homepage = "https://pythonhosted.org/pyDOE/";
      license = licenses.mit;
      maintainers = with lib.maintainers; [ malcolms ];
  };
}

