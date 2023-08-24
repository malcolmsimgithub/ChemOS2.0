{ buildPythonPackage
, fetchPypi
, pytest
, lib
, numpy
}:

buildPythonPackage rec {
  pname = "matter-chimera";
  version = "1.0";
  format = "setuptools";

  src = fetchPypi {
    inherit pname version format;
    sha256 = "9f20eb02c4997cdcea551d21c88dbc00b82cba28cdaf76ef6fa6599cec7720f9";
  };

  propagatedBuildInputs = [
    numpy
  ];

  meta = with lib; {
    description = "Chimera is a general purpose achievement scalarizing function for multi-objective optimization.";
    homepage = "https://github.com/aspuru-guzik-group/chimera";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

