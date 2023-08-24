{ buildPythonPackage
, black
, fetchurl
, flake8
, gurobipy
, isort
, lib
, numpy
, nbsphinx
, nbval
, opt-einsum
, pandas
, pypandoc
, pyro-api
, pytest
, seaborn
, scikit-learn
, scipy
, sphinx
, sqlalchemy
, torch
, tqdm
, yapf
}:
  
buildPythonPackage rec {
  pname = "pyro-ppl";
  version = "1.8.4";
  format = "wheel";

  src = fetchurl{
      url = "https://files.pythonhosted.org/packages/b5/b1/ccceeae368b7e2b5504229e74ad584e4b8071faeef23b0e888d1c9d8ef3d/pyro_ppl-1.8.4-py3-none-any.whl";
      sha256 = "294f78f28f2fe7bbea2792bd6bd8c69b7cfe493cf8940cac97a9b5d0e7f194cd";
  };
  doCheck = true;

  propagatedBuildInputs = [
    torch
    numpy
    pyro-api
    tqdm
    opt-einsum
    black
    flake8
    isort
    nbsphinx
    nbval
    pypandoc
    pytest
    scipy
    sphinx
    yapf
  ];

  meta = with lib; {
    description = "Deep Universal Probabilistic Programming";
    homepage = "https://pyro.ai";
    license = licenses.asl20;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

