{ pkgs ? import <nixpkgs> { } }:

with pkgs;

let

  singleton-decorator = python310Packages.buildPythonPackage rec {
    pname = "singleton-decorator";
    version = "1.0.0";
    format = "setuptools";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/33/98/a8b5c919bee1152a9a1afd82014431f8db5882699754de50d1b3aba4d136/singleton-decorator-1.0.0.tar.gz";
      sha256 = "1a90ad8a8a738be591c9c167fdd677c5d4a43d1bc6b1c128227be1c5e03bee07";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
    ];
  };

  pyunpack = python310Packages.buildPythonPackage rec {
    pname = "pyunpack";
    version = "0.3";
    format = "wheel";
    src = fetchurl {
      url ="https://files.pythonhosted.org/packages/72/20/5a6dcb0d28529ce6efe850755994c80817279eecf08620003775fda3b914/pyunpack-0.3-py2.py3-none-any.whl";
      sha256 ="8f517cfc71215f37f74cf3a7668028828c68dc76f4d02e7a69f227ce978d51a3";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
      entrypoint2
      EasyProcess
    ];
  };
  blinker = python310Packages.buildPythonPackage rec {
    pname = "blinker";
    version = "1.6.2";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/0d/f1/5f39e771cd730d347539bb74c6d496737b9d5f0a53bc9fdbf3e170f1ee48/blinker-1.6.2-py3-none-any.whl";
      sha256 = "c3d739772abb7bc2860abf5f2ec284223d9ad5c76da018234f6f50d6f31ab1f0";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
    ];
  };
  cachetools = python310Packages.buildPythonPackage rec {
    pname = "cachetools";
    version = "5.3.0";
    format = "wheel";
    src = fetchurl {
      url= "https://files.pythonhosted.org/packages/db/14/2b48a834d349eee94677e8702ea2ef98b7c674b090153ea8d3f6a788584e/cachetools-5.3.0-py3-none-any.whl";
      sha256 = "429e1a1e845c008ea6c85aa35d4b98b65d6a9763eeef3e37e92728a12d1de9d4";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
    ];
  };
  protobuf = python310Packages.buildPythonPackage rec {
    pname = "protobuf";
    version = "3.20.3";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/8d/14/619e24a4c70df2901e1f4dbc50a6291eb63a759172558df326347dce1f0d/protobuf-3.20.3-py2.py3-none-any.whl";
      sha256 = "a7ca6d488aa8ff7f329d4c545b2dbad8ac31464f1d8b1c87ad1346717731e4db";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
    ];
  };
  pympler = python310Packages.buildPythonPackage rec {
    pname = "pympler";
    version = "1.0.1";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/2c/42/41e1469ed0b37b9c8532cb8074bea179f7d85ee7e82a59b5b6c289ed6045/Pympler-1.0.1-py3-none-any.whl";
      sha256 = "d260dda9ae781e1eab6ea15bacb84015849833ba5555f141d2d9b7b7473b307d";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
    ];
  };
  streamlit = python310Packages.buildPythonPackage rec {
    pname = "streamlit";
    version = "1.22.0";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/b1/26/2add66d2e2febd6b05efd48ce02fc2d979d805b844143f2fa9fd7e867ade/streamlit-1.22.0-py2.py3-none-any.whl";
      sha256 = "520dd9b9e6efb559b5a9a22feadb48b1e6f0340ec83da3514810059fdecd4167";
    };

    doCheck = false;

    propagatedBuildInputs = with python310Packages; [
      packaging
      pympler
      protobuf
      cachetools
      blinker
      requests
      pandas
      toml
      pyarrow
      pillow
      validators
      rich
      click
      importlib-metadata
      tenacity
      pydeck
      GitPython
      tzlocal
      tornado
      watchdog
      pyunpack
      altair
    ];
  };
  sila2 = python310Packages.buildPythonPackage rec {
    pname = "sila2";
    version = "0.10.1";
    format = "wheel";
    src = fetchurl {
      url = "https://files.pythonhosted.org/packages/66/73/a1bcc3ed7e27e06531623751d4ef40bfed0848aaab8018be82a3efdb9345/sila2-0.10.1-py3-none-any.whl";
      sha256 = "2bae6eacaa5ea47d0f625bf81c50332b5a733c3af3703033084cf1441fde99f0";
    };

    doCheck = true;

    propagatedBuildInputs = with python310Packages; [
      lxml
      typing-extensions
      cryptography
      grpcio-tools
      zeroconf
      black
      flake8
      pytest
      pytest-cov
      sqlalchemy
    ];
  pythonImportsCheck = [ "sila2" ];

  meta = with lib; {
    description = "SILA2 bindings for python";
    homepage = "https://sila-standard.com";
    license = with licenses; [ mit ];
};
};
  linear_operator = python310Packages.buildPythonPackage rec {
      pname = "linear_operator";
      version = "0.3.0";
      format = "setuptools";

      src = python310Packages.fetchPypi {
        inherit pname version format;
        

        sha256 = "sha256-hL9XJjGn4Vdt5pINgWAMoP7c9r2i8p269EDW5yzmq6s="; # TODO
      };
      
    propagatedBuildInputs = with python310Packages; [
      torch
      scipy
    ];
    };      
  gpytorch = python310Packages.buildPythonPackage rec {
      pname = "gpytorch";
      version = "1.9.1";
      format = "setuptools";

      src = python310Packages.fetchPypi {
        inherit pname version format;
        

        sha256 = "sha256-C9u6b21ZV6D0Pvbcf+w5xH6KVfYyyjN2DGGJ8lmzzMM="; # TODO
      };
      
    propagatedBuildInputs = with python310Packages; [
      torch
      scikit-learn
      linear_operator
    ];
    };      
    
  pyro-ppl= python310Packages.buildPythonPackage rec {
      pname = "pyro-ppl";
      version = "1.8.4";
      format = "wheel";

      src = fetchurl{
          url = "https://files.pythonhosted.org/packages/b5/b1/ccceeae368b7e2b5504229e74ad584e4b8071faeef23b0e888d1c9d8ef3d/pyro_ppl-1.8.4-py3-none-any.whl";
          sha256 = "294f78f28f2fe7bbea2792bd6bd8c69b7cfe493cf8940cac97a9b5d0e7f194cd";
      };
      doCheck = true;
    propagatedBuildInputs = with python310Packages; [
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
    };      

    botorch = python310Packages.buildPythonPackage rec {
      pname = "botorch";
      version = "0.8.1";
      format = "wheel";

      src = fetchurl {
        url = "https://files.pythonhosted.org/packages/36/b5/679cff67119aa08f45fc7b25318424dcbd3f9c207a671212801d02bf3fd3/botorch-0.8.1-py3-none-any.whl";
        sha256 = "f1f89a73c51d5e5db5e258443d8125e5db535dd8ce0f32e11513bc7f1cb28fe7";
      };
      
    propagatedBuildInputs = with python310Packages; [
      torch
      scikit-learn
      linear_operator
      multipledispatch
      scipy
      gpytorch
      pyro-ppl
    ];
    };      

    matter-golem = python310Packages.buildPythonPackage rec {
      pname = "matter-golem";
      version = "1.0";
      format = "setuptools";

      src = python310Packages.fetchPypi {
        inherit pname version format;
        sha256 = "9df35debdc40d70f724d644bfa755654e63c33279690a28c2d71e33294ee16b4";
      };
      
    propagatedBuildInputs = with python310Packages; [
      pandas
      numpy
      scipy
      scikit-learn
      deap
      pytest
    ];
    };      

    matter-chimera = python310Packages.buildPythonPackage rec {
      pname = "matter-chimera";
      version = "1.0";
      format = "setuptools";

      src = python310Packages.fetchPypi {
        inherit pname version format;
        sha256 = "9f20eb02c4997cdcea551d21c88dbc00b82cba28cdaf76ef6fa6599cec7720f9"; # TODO
      };
      
    propagatedBuildInputs = with python310Packages; [
      numpy
    ];
  };

  atlas = python310Packages.buildPythonPackage rec {
        name = "atlas";
        format = "setuptools";
        src = fetchGit {
          name = "atlas";
          url = "file:///home/malcolm/env/atlas-dynamic_search";
          rev = "0b944a973ffdb6002ccfc998d81695cd2b1b868a";
        };

        doCheck = false;

        propagatedBuildInputs = with python310Packages; [
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
          gspread
        ];
  };

  gurobipy = python310Packages.buildPythonPackage rec{
      pname = "gurobipy";
      format = "wheel";
      version = "10.0.1";
      src = fetchurl {

        url = "https://files.pythonhosted.org/packages/0f/42/4e007957b38fcd1dfe745c32cef23b51dd4a68719544ddb657d18b816724/gurobipy-10.0.1-cp310-cp310-manylinux2014_x86_64.whl";
        sha256 = "ac0da9ae09c681ef59528b5eff28a60c858825c5976dac2fa6a1f5c0b7a324ad";

      };


    };

  olympus = python310Packages.buildPythonPackage rec{
        pname = "olympus";
        format = "setuptools";
        version = "0.0.1b0";
        src = fetchGit {
            url = "file:///home/malcolm/env/olympus";
            rev = "8945b0b8067a014359c095942754d1f923ca9075";
          };

        propagatedBuildInputs = with python310Packages; [
          pandas
          numpy
          pytest
          sqlalchemy
          seaborn
          gurobipy
          scikit-learn
        ];


  };

  notification-apps = python310Packages.buildPythonPackage rec{
        pname = "notification-apps";
        format = "setuptools";
        version = "0.0.1b0";
        src = fetchGit {
            url = "file:///home/malcolm/notification-apps";
            rev = "1c681922fb9556ec34e01aeaad1452dd2d62e614";   
          };

        propagatedBuildInputs = with python310Packages; [
        slack-sdk
        ];


  };
  molar = python310Packages.buildPythonPackage rec {
      pname = "molar";
      format = "wheel";
      version = "0.4.4";
      src = fetchurl {
  
        url = "https://files.pythonhosted.org/packages/98/38/6ea4704b1ce80b2e0f70a10fb21f7e46c5b31c18cca245430b1581c13eda/molar-0.4.4-py3-none-any.whl";
        sha256 = "ec70d75c94ba1c3b499a719924c329bd18806d842041860f52fcec86d43ebea9";
      };

      propagatedBuildInputs = with python310Packages; [
        pandas
        numpy
        rich
        pydantic
        python-jose
        fastapi
        magic
        click

    ];
  };
in
  mkShell { buildInputs = with python310Packages; [ notification-apps sila2 packaging streamlit singleton-decorator python-telegram-bot tkinter pandas typer gpytorch botorch matter-chimera matter-golem olympus atlas molar psycopg2 rdkit pyunpack patool];

  shellHook = ''
    export MOLAR_PASSWORD=sefvac-8mejvu-wujziJ
    export MOLAR_USER=malcolm.sim@mail.utoronto.ca
    export OLYMPUS_SCRATCH="/home/malcolm/.scratch"
    export SLACK_API_TOKEN="xoxb-2544853610-5246320937812-8EmeFWdEkMS9Q0YaZrZjRcSO"
  '';
}
