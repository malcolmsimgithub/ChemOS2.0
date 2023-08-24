#with import <nixpkgs> {};
#with pkgs.python310Packages;
{ buildPythonPackage
, gpytorch
, multipledispatch
, lib
, linear_operator
, fetchurl
, scikit-learn
, scipy
, pyro-ppl
, torch
}:

buildPythonPackage rec {
  pname = "botorch";
  version = "0.8.1";
  format = "wheel";

  src = fetchurl {
    url = "https://files.pythonhosted.org/packages/36/b5/679cff67119aa08f45fc7b25318424dcbd3f9c207a671212801d02bf3fd3/botorch-0.8.1-py3-none-any.whl";
    sha256 = "f1f89a73c51d5e5db5e258443d8125e5db535dd8ce0f32e11513bc7f1cb28fe7";
  };

  propagatedBuildInputs = [
    torch
    scikit-learn
    linear_operator
    multipledispatch
    scipy
    gpytorch
    pyro-ppl
  ];

  meta = with lib; {
    description = "BoTorch is a library for Bayesian Optimization built on PyTorch.";
    homepage = "https://botorch.org";
    license = licenses.mit;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

