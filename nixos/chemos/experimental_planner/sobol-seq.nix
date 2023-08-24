{ buildPythonPackage
, fetchurl
, lib
, numpy
, scipy
}:

buildPythonPackage rec {
  pname = "sobol-seq";
  version = "0.2.0";
  format = "wheel";

  doCheck=false;

  src = fetchurl {
    url = "https://files.pythonhosted.org/packages/e4/df/6c4ad25c0b48545a537b631030f7de7e4abb939e6d2964ac2169d4379c85/sobol_seq-0.2.0-py3-none-any.whl";
    sha256 = "277ab767250a20b440fc74df8b6f4d79773949d5770927e1cee83e8de026b704";
  };

  propagatedBuildInputs = [
    numpy
    scipy
  ];

  meta = with lib; {
      description = "Sobol sequence generator";
      homepage = "https://pypi.org/project/sobol-seq/";
      license = licenses.mit;
      maintainers = with lib.maintainers; [ malcolms ];
  };
}

