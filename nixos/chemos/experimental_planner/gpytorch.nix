{ buildPythonPackage
, lib
, linear_operator
, fetchPypi
, fetchurl
, torch
, scikit-learn
}:

buildPythonPackage rec {
  pname = "gpytorch";
  version = "1.9.1";
  format = "wheel";
  doCheck = true;
  src = fetchurl {
    url = "https://files.pythonhosted.org/packages/0b/8d/2baa077efc7ad94fa47b8a80324166538f75d76cb4ac2c7833bbc2938cc9/gpytorch-1.9.1-py3-none-any.whl";
    sha256 = "sha256-/u7bn9QjCopr9+0y/3SaDef75S7XIxigzizIXz2dSag=";
    };

  propagatedBuildInputs = [
    torch
    scikit-learn
    linear_operator
  ];

  meta = with lib; {
    description = "gpytorch";
    homepage = "https://gpytorch.ai";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

