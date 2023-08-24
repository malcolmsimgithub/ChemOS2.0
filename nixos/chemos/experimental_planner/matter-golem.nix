{ buildPythonPackage
, deap
, fetchPypi
, pytest
, lib
, pandas
, numpy
, scikit-learn
, scipy
}:
  
buildPythonPackage rec {
  pname = "matter-golem";
  version = "1.0";
  format = "setuptools";

  src = fetchPypi {
    inherit pname version format;
    sha256 = "9df35debdc40d70f724d644bfa755654e63c33279690a28c2d71e33294ee16b4";
  };

   propagatedBuildInputs = [
     pandas
     numpy
     scipy
     scikit-learn
     deap
     pytest
  ];
  
  meta = with lib; {
    description = "Golem: An Algorithm for Robust Experiment and Process Optimization";
    homepage = "https://pypi.org/project/matter-golem/";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

