{ buildPythonPackage
, lib
, fetchPypi
, torch
, scipy
}:

buildPythonPackage rec {
  pname = "linear_operator";
  version = "0.3.0";
  format = "setuptools";

  src = fetchPypi {
    inherit pname version format;
    sha256 = "sha256-hL9XJjGn4Vdt5pINgWAMoP7c9r2i8p269EDW5yzmq6s=";
  };

  propagatedBuildInputs = [
    torch
    scipy
  ];

  meta = with lib; {
    description = "A linear operator implementation, primarily designed for finite-dimensional positive definite operators (i.e. kernel matrices).";
    homepage = "https://pypi.org/project/linear-operator/";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}
