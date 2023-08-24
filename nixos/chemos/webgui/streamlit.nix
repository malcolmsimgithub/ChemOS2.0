{ buildPythonPackage
, lib
, fetchurl
, packaging
, pympler
, protobuf
, cachetools
, blinker
, requests
, pandas
, toml
, pyarrow
, pillow
, validators
, rich
, click
, importlib-metadata
, tenacity
, pydeck
, GitPython
, tzlocal
, tornado
, watchdog
, altair
}:
  
buildPythonPackage rec {
  pname = "streamlit";
  version = "1.24.1";
  format = "wheel";
  src = fetchurl {
    url = "https://files.pythonhosted.org/packages/6a/23/a788a838a3466ac009012a998e73e7cb3eab896a26b4eab7bcb2195ef9e0/streamlit-1.24.1-py2.py3-none-any.whl";
    sha256 = "569211b426ca078c0c2959a6c9a1413c2e81eca23e769fbf12308ba5bffd1f49";
  };

  doCheck = false;

  propagatedBuildInputs = [
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
    altair
  ];


  meta = with lib; {
    description = "A faster way to build and share data apps";
    homepage = "https://streamlit.io";
    license = licenses.asl20;
    maintainers = with lib.maintainers; [ malcolms ];
  };
}

