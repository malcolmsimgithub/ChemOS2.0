{ lib
, buildPythonPackage
, fetchurl
, pytest
, pandas
, numpy
, scipy
, seaborn
, gurobipy
, scikit-learn
, sqlalchemy
}:                                                           
                                                               
buildPythonPackage rec {
    pname = "singleton-decorator";
    version = "1.0.0";
    format = "setuptools";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/33/98/a8b5c919bee1152a9a1afd82014431f8db5882699754de50d1b3aba4d136/singleton-decorator-1.0.0.tar.gz";
      sha256 = "1a90ad8a8a738be591c9c167fdd677c5d4a43d1bc6b1c128227be1c5e03bee07";
    };

    doCheck = true;

    propagatedBuildInputs = [
      pandas
      numpy
      pytest
      sqlalchemy
      seaborn
      gurobipy
      scikit-learn
    ];

    pythonImportsCheck = ["singleton_decorator"];

    meta = with lib; {
      description = "Python singleton decorator";
      homepage = "https://pypi.org/project/singleton-decorator/";
      license = licenses.gpl3;
      maintainers = with lib.maintainers; [ malcolms ];
    };
}
