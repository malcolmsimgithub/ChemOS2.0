#with import <nixpkgs> {};
#with pkgs.python310Packages;
{ buildPythonPackage
, lxml
, typing-extensions
, cryptography
, grpcio-tools
, zeroconf
, black
, flake8
, pytest
, pytest-cov
, lib
, pythonOlder
, fetchurl
, periodictable
, pandas
, numpy
, scipy
}:

buildPythonPackage rec {
    pname = "sila2";
    version = "0.10.1";
    format = "wheel";

    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/66/73/a1bcc3ed7e27e06531623751d4ef40bfed0848aaab8018be82a3efdb9345/sila2-0.10.1-py3-none-any.whl";
      sha256 = "2bae6eacaa5ea47d0f625bf81c50332b5a733c3af3703033084cf1441fde99f0";
    };

    doCheck = true;

    propagatedBuildInputs = [
      lxml
      typing-extensions
      cryptography
      grpcio-tools
      zeroconf
      black
      flake8
      pytest
      pytest-cov
      periodictable
      scipy
      numpy
      pandas
    ];

  pythonImportsCheck = [ "sila2" ];

  meta = with lib; {
    description = "SILA2 bindings for python";
    homepage = "https://sila-standard.com";
    license = licenses.mit;
  };
}
