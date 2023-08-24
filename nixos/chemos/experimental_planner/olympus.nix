{ lib
, buildPythonPackage
, pytest
, fetchFromGitHub
, pandas
, numpy
, scipy
, gurobipy
, seaborn
, scikit-learn
, sqlalchemy
}:
                                                               
buildPythonPackage rec {                                      
  pname = "olympus";
  format = "setuptools";
  version = "0.0.1b0";
  src = fetchFromGitHub {
    owner = "aspuru-guzik-group";
    repo = "olympus";
    rev = "440b6b58ebfcaa2391cff7e94b570fb4fda98d68";
    sha256 = "sha256-+uyy4XQIKrHwIsKU/n1tOllr0J9A/N1kJz/oso5PJXc=";
  };

  preConfigure = ''
    sed -i 's/"SQLAlchemy==1.4.45"/"SQLAlchemy"/' setup.py
  '';

  propagatedBuildInputs = [
    pandas
    numpy
    pytest
    sqlalchemy
    seaborn
    gurobipy
    scikit-learn
  ];

  pythonImportsCheck = ["olympus"];

  meta = with lib; {
    description = "A Benchmarking framework for Bayesian optimization in chemsitry";
    homepage = "https://github.com/aspuru-guzik-group/olympus";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

