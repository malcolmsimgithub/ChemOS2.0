{ lib
, buildPythonPackage
, python310Packages
, fetchurl
, EasyProcess
, entrypoint2
}:
                                                               
buildPythonPackage rec {
  pname = "pyunpack";
  version = "0.3";
  format = "wheel";
  src = fetchurl {
    url ="https://files.pythonhosted.org/packages/72/20/5a6dcb0d28529ce6efe850755994c80817279eecf08620003775fda3b914/pyunpack-0.3-py2.py3-none-any.whl";
    sha256 ="8f517cfc71215f37f74cf3a7668028828c68dc76f4d02e7a69f227ce978d51a3";
  };

  doCheck = true;

  propagatedBuildInputs = with python310Packages; [
      EasyProcess
      entrypoint2
  ];

  pythonImportsCheck = ["pyunpack"];

  meta = with lib; {
    description = "unpack archive files in Python";
    homepage = "https://github.com/ponty/pyunpack";
    license = licenses.bsd2;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}
