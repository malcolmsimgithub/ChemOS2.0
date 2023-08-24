{ buildPythonPackage
, botorch
, deap
, fetchFromGitHub
, gspread
, ipykernel
, jupyterlab
, matplotlib
, matter-chimera
, matter-golem
, pyDOE
, pymoo
, python310Packages
, pytest
, lib
, pandas
, rich
, scipy
, seaborn
, scikit-learn
, sobol-seq
, sqlalchemy
, tqdm
}:

buildPythonPackage rec {
  name = "atlas";
  format = "setuptools";
  src = fetchFromGitHub {
    owner = "rileyhickman";
    repo = "atlas";
    rev = "9aa18933bcfc9e1d363dfc9e15e0cbdcc0aa7e01";
    sha256 = "sha256-mZ39lJoP1EoFBxC+upJhGGWFhKXgep5jQUDf0hExWN0=";
  };

  doCheck = false;

  propagatedBuildInputs = [
    botorch
    matter-chimera
    matter-golem
    scikit-learn
    scipy
    tqdm
    ipykernel
    jupyterlab
    matplotlib
    seaborn
    pandas
    sqlalchemy
    deap
    rich
    pytest
    pymoo
    gspread
    pyDOE
    sobol-seq
  ];
  meta = with lib; {
      description = "A brain for self-driving laboratories";
      homepage = "https://github.com/aspuru-guzik-group/atlas";
      license = licenses.mit;
      maintainers = with lib.maintainers; [ malcolms ];
  };
}

